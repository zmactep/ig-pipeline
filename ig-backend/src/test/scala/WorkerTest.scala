import akka.actor._
import akka.testkit.{TestKit, ImplicitSender}
import akka.util.Timeout
import akka.pattern.{pipe, ask}
import scala.concurrent.duration._
import master.Master
import org.scalatest.{WordSpec, BeforeAndAfterAll}
import org.scalatest.matchers.MustMatchers
import scala.concurrent.{ExecutionContext, Await, Future}
import workers.{SimpleWorker, Worker}

class WorkerTest extends TestKit(ActorSystem("WorkerTest"))
with ImplicitSender
with WordSpec
with BeforeAndAfterAll
with MustMatchers {

  implicit val askTimeout = Timeout(5 second)
  override def afterAll() {
    system.shutdown()
  }

  val master = system.actorOf(Props[Master], "master")
  val w1 = worker("master")
  val w2 = worker("master")
  val w3 = worker("master")

  "Worker" should {
    "work" in {
      master ! "{\"task\" : \"model list\", \"input\": {\"group\": \"regions\"}}"
      master ! "{\"task\" : \"model list\", \"input\": {\"group\": \"regions\"}}"
      expectMsgAllOf("Your job is enqueued with id = 0", "Your job is enqueued with id = 1")

      Thread.sleep(2000)

      master ! "{\"result_for\":\"0\"}"
      master ! "{\"result_for\":\"1\"}"
      master ! "{\"result_for\":\"2\"}"
      expectMsgAllOf("{\"status\": \"ok\"}", "{\"status\": \"ok\"}", "Your job id not found")
    }

    "work with Futures" in {
      import ExecutionContext.Implicits.global
      val task = "{\"task\" : \"model list\", \"input\": {\"group\": \"regions\"}}"
      val fs1 = Future.sequence(List(task, task).map { s => master ? s })
      Await.result(fs1, 1 second) must be (List("Your job is enqueued with id = 2", "Your job is enqueued with id = 3"))

      Thread.sleep(2000)

      val fs2 = Future.sequence(List("{\"result_for\":\"2\"}", "{\"result_for\":\"3\"}", "{\"result_for\":\"42\"}").map { s => master ? s })
      Await.result(fs2, 1 second) must be (List("{\"status\": \"ok\"}", "{\"status\": \"ok\"}", "Your job id not found"))
    }
  }
  private def worker(name: String) = system.actorOf(Props(
    new SimpleWorker(ActorPath.fromString(
      "akka://%s/user/%s".format(system.name, name)))))
}
