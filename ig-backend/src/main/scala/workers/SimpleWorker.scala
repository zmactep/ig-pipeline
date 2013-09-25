package workers

import sys.process._
import scala.concurrent.Future
import akka.actor.{ActorPath, ActorRef, ActorLogging}
import akka.pattern.pipe

/**
 * Created with IntelliJ IDEA.
 * User: Kos
 * Date: 25.09.13
 * Time: 10:02
 * To change this template use File | Settings | File Templates.
 */

class SimpleWorker(masterLocation: ActorPath) extends Worker(masterLocation) with ActorLogging{
  // We'll use the current dispatcher for the execution context.
  // You can use whatever you want.
  implicit val ec = context.dispatcher

  def doWork(workSender: ActorRef, msg: Any): Unit = {
    Future {
      val cwd = context.system.settings.config.getString("ig-backend.tools_root")
      val execResult = Process(Seq("./train_model_example.sh"), new java.io.File(cwd)).!!
      workSender ! "Your request: " + msg + ". Response: " + execResult
      WorkComplete("done")
    } pipeTo self
  }
}

