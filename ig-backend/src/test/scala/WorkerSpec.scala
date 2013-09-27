import akka.actor.{ActorSystem, ActorPath, ActorRef, Props}
import akka.testkit.{TestKit, ImplicitSender}
import akka.util.Timeout
import akka.pattern.{pipe, ask}
import scala.concurrent.duration._
import master.Master
import org.scalatest.{WordSpec, BeforeAndAfterAll}
import org.scalatest.matchers.MustMatchers
import scala.concurrent.{ExecutionContext, Await, Future}
import workers.Worker

class TestWorker(masterLocation: ActorPath) extends Worker(masterLocation) {
  // We'll use the current dispatcher for the execution context.
  // You can use whatever you want.
  implicit val ec = context.dispatcher

  def doWork(workSender: ActorRef, msg: Any): Unit = {
    Future {
      workSender ! msg
      WorkComplete("done")
    } pipeTo self
  }
}

class BadTestWorker(masterLocation: ActorPath) extends Worker(masterLocation) {
  // We'll use the current dispatcher for the execution context.
  // You can use whatever you want.
  implicit val ec = context.dispatcher

  def doWork(workSender: ActorRef, msg: Any): Unit = context.stop(self)
}

class WorkerSpec extends TestKit(ActorSystem("WorkerSpec"))
with ImplicitSender
with WordSpec
with BeforeAndAfterAll
with MustMatchers {

  implicit val askTimeout = Timeout(1 second)
  override def afterAll() {
    system.shutdown()
  }

  def worker(name: String) = system.actorOf(Props(
    new TestWorker(ActorPath.fromString(
      "akka://%s/user/%s".format(system.name, name)))))

  def badWorker(name: String) = system.actorOf(Props(
    new BadTestWorker(ActorPath.fromString(
      "akka://%s/user/%s".format(system.name, name)))))

  "Worker" should {
    "work" in {
      // Spin up the master
      val m = system.actorOf(Props[Master], "master")
      // Create three workers
      val w1 = worker("master")
      val w2 = worker("master")
      val w3 = worker("master")
      // Send some work to the master
      m ! "Hithere"
      m ! "Guys"
      // We should get it all back
      expectMsgAllOf("Your job is enqueued with id = 0", "Your job is enqueued with id = 1")
      m ! "get_result: 0"
      m ! "get_result: 1"
      expectMsgAllOf("Your result: done", "Your result: done")
    }
    "still work if one dies" in { //{2
    // Spin up the master
    val m = system.actorOf(Props[Master], "master2")
      // Create three workers
      val w1 = worker("master2")
      val w2 = badWorker("master2")
      // Send some work to the master
      m ! "Hithere"
      m ! "Guys"
      // We should get it all back
      expectMsgAllOf("Your job is enqueued with id = 0", "Your job is enqueued with id = 1")
      m ! "get_result: 0"
      m ! "get_result: 2"
      //TODO fix a bug with task_id reassignment
      expectMsgAllOf("Your result: done", "Your job id not found")
    } //}2
    "work with Futures" in { //{2
      import ExecutionContext.Implicits.global
    // Spin up the master
    val m = system.actorOf(Props[Master], "master3")
      // Create three workers
      val w1 = worker("master3")
      val w2 = worker("master3")
      val w3 = worker("master3")
      val fs1 = Future.sequence(List("Hithere", "Guys").map { s => m ? s })
      Await.result(fs1, 1 second) must be (List("Your job is enqueued with id = 0", "Your job is enqueued with id = 1"))
      // We should get it all back
      val fs2 = Future.sequence(List("get_result: 0", "get_result: 1").map { s => m ? s })
      Await.result(fs2, 1 second) must be (List("Your result: done", "Your result: done"))
    } //}2
  }
}
