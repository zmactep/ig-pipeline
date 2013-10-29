import akka.actor._
import akka.testkit.{TestKit, ImplicitSender}
import akka.util.Timeout
import akka.pattern.ask
import org.scalatest.mock.MockitoSugar
import scala.concurrent.duration._
import master.Master
import org.scalatest.{WordSpec, BeforeAndAfterAll}
import org.scalatest.matchers.MustMatchers
import scala.concurrent.{ExecutionContext, Await, Future}
import workers.{SimpleWorker}

import org.mockito.Mockito._
import org.mockito.Matchers._
import java.sql._


class WorkerTest extends TestKit(ActorSystem("WorkerTest"))
with ImplicitSender
with MockitoSugar
with WordSpec
with BeforeAndAfterAll
with MustMatchers {

  implicit val askTimeout = Timeout(5 second)
  override def afterAll() {
    system.shutdown()
  }

  val connection = mock[Connection]
  val statement = mock[Statement]
  val preparedStatement = mock[PreparedStatement]
  val resultSet = mock[ResultSet]
  when(connection.createStatement()).thenReturn(statement)
  when(connection.prepareStatement(anyString)).thenReturn(preparedStatement)

  val master = system.actorOf(Master.props(connection), "master")
  val w1 = worker("master")
  val w2 = worker("master")
  val w3 = worker("master")

  when(resultSet.next()).thenReturn(true).thenReturn(true).thenReturn(true).thenReturn(true).
    thenReturn(false).thenReturn(true).thenReturn(true).thenReturn(true).thenReturn(true).thenReturn(false)
  when(resultSet.getInt(1)).thenReturn(0).thenReturn(1).thenReturn(2).thenReturn(3)
  when(resultSet.getString("status")).thenReturn("ok").thenReturn("ok")
  when(resultSet.getString("result")).thenReturn("result1").thenReturn("result2").thenReturn("result3").thenReturn("result4")
  when(preparedStatement.executeUpdate()).thenReturn(1)
  when(preparedStatement.executeQuery()).thenReturn(resultSet)
  when(statement.executeQuery(anyString)).thenReturn(resultSet)

  "Worker" should {
    "work" in {

      master ! "{\"task\" : \"test\", \"input\": {\"group\": \"regions\"}}"
      master ! "{\"task\" : \"test\", \"input\": {\"group\": \"regions\"}}"
      expectMsgAllOf("{\"id\":0}", "{\"id\":1}")

      Thread.sleep(2000)

      master ! "{\"result_for\":\"0\"}"
      master ! "{\"result_for\":\"1\"}"
      master ! "{\"result_for\":\"2\"}"
      expectMsgAllOf("{\"status\": \"ok\", data: [result1]}", "{\"status\": \"ok\", data: [result2]}", "Your job id not found")
    }

    "work with Futures" in {
      import ExecutionContext.Implicits.global

      val task = "{\"task\" : \"test\", \"input\": {\"group\": \"regions\"}}"
      val fs1 = Future.sequence(List(task, task).map { s => master ? s })
      Await.result(fs1, 1 second) must be (List("{\"id\":2}", "{\"id\":3}"))

      Thread.sleep(2000)

      val fs2 = Future.sequence(List("{\"result_for\":\"2\"}", "{\"result_for\":\"3\"}", "{\"result_for\":\"42\"}").map { s => master ? s })
      Await.result(fs2, 1 second) must be (List("{\"status\": \"ok\", data: [result3]}", "{\"status\": \"ok\", data: [result4]}", "Your job id not found"))
    }
  }
  private def worker(name: String) = system.actorOf(Props(
    new SimpleWorker(ActorPath.fromString(
      "akka://%s/user/%s".format(system.name, name)))))
}
