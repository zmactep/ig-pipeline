package master

import akka.actor._
import org.json.JSONObject
import scala.Tuple2
import akka.actor.Terminated
import scala.Some
import utils.DbUtils
import java.sql.{Connection}

/**
 * Created with IntelliJ IDEA.
 * User: Kos
 * Date: 25.09.13
 * Time: 11:35
 * To change this template use File | Settings | File Templates.
 */

object Master {
  def props(connection: Connection): Props = Props(new Master(connection))
}

//https://github.com/derekwyatt/akka-worker-pull
class Master(conn: Connection) extends Actor with ActorLogging {
  import MasterWorkerProtocol._
  import scala.collection.mutable.{Map, Queue}

  classOf[com.mysql.jdbc.Driver]
  val prepInsert = conn.prepareStatement("INSERT INTO ig_backend_tasks (params, result, status) VALUES (?, ?, ?)")
  val prepUpdate = conn.prepareStatement("UPDATE ig_backend_tasks set result = ?, status = ? WHERE id = ?")
  val prepSelect = conn.prepareStatement("SELECT result, status from  ig_backend_tasks WHERE id = ?")

  // Holds known workers and what they may be working on
  //jobId, worksender, work
  val workers = Map.empty[ActorRef, Option[Tuple2[Int, Tuple2[ActorRef, Any]]]]

  // Holds the incoming list of work to be done as well
  // as the memory of who asked for it
  val workQ = Queue.empty[Tuple2[ActorRef, Any]]

  // Notifies workers that there's work available, provided they're
  // not already working on something
  def notifyWorkers(): Unit = {
    if (!workQ.isEmpty) {
      workers.foreach {
        case (worker, m) if (m.isEmpty) => worker ! WorkIsReady
        case _ =>
      }
    }
  }

  def receive = {
    // Worker is alive. Add him to the list, watch him for
    // death, and let him know if there's work to be done
    case WorkerCreated(worker) =>
      log.info("Worker created: {}", worker)
      context.watch(worker)
      workers += (worker -> None)
      notifyWorkers()

    // A worker wants more work.  If we know about him, he's not
    // currently doing anything, and we've got something to do,
    // give it to him.
    case WorkerRequestsWork(worker) =>
      log.info("Worker requests work: {}", worker)
      if (workers.contains(worker)) {
        if (workQ.isEmpty)
          worker ! NoWorkToBeDone
        else if (workers(worker) == None) {
          val (workSender, work) = workQ.dequeue()
          val jobId = DbUtils.saveWorkTask(work, prepInsert, conn)
          workers += (worker -> Some(jobId -> (workSender -> work)))
          worker ! WorkToBeDone(work, jobId)
          workSender ! new JSONObject().put("id", jobId).toString
        }
      }

    // Worker has completed its work and we can clear it out
    case WorkIsDone(worker, result) =>
      if (!workers.contains(worker))
        log.error("Blurgh! {} said it's done work but we didn't know about him", worker)
      else {
        val (jobId, _) = workers(worker).get
        workers += (worker -> None)
        DbUtils.updateTask(jobId, result, prepUpdate, conn)
        log.debug("Your job result #{} is saved", jobId)
      }

    // A worker died.  If he was doing anything then we need
    // to give it to someone else so we just add it back to the
    // master and let things progress as usual
    case Terminated(worker) =>
      if (workers.contains(worker) && workers(worker) != None) {
        log.error("Blurgh! {} died while processing {}", worker, workers(worker))
        // Send the work that it was doing back to ourselves for processing
        val (_, (workSender, work)) = workers(worker).get   //TODO jobId is changed here - we will not find task by old id
        self.tell(work, workSender)
      }
      workers -= worker

    // Anything other than our own protocol is "work to be done"
    case work =>
      val cmd = work.toString.trim
      if (cmd.contains("result_for")) {
        try {
          val jobId = new JSONObject(cmd).getString("result_for").toInt
          log.info("Requesting result for {}", jobId.toInt)
          sender ! DbUtils.findTask(jobId, prepSelect, conn)
        } catch {
          case e: Exception => sender ! "Missing job ID in query"
        }

      } else {
        log.info("Queueing {}", work)
        workQ.enqueue(sender -> work)
        notifyWorkers()
      }
  }


}
