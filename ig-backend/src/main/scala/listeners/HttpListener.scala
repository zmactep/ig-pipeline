package listeners

import akka.actor.{ Actor, ActorLogging, ActorRef, Props, Terminated }
import akka.io.{ IO, Tcp }
import akka.pattern.ask
import akka.util.Timeout
import java.net.InetSocketAddress
import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Await}
import spray.can.Http
import spray.http.HttpMethods._
import spray.http._
import spray.httpx.unmarshalling._
import spray.can.server.ServerSettings

/**
 * Created with IntelliJ IDEA.
 * User: Kos
 * Date: 25.09.13
 * Time: 8:34
 * To change this template use File | Settings | File Templates.
 */

object HttpListener {

  def props(host: String, port: Int): Props =
    Props(new HttpListener(host, port))
}

class HttpListener(host: String, port: Int) extends Actor with ActorLogging {

  import context.system

  IO(Http) ! Http.Bind(self, host, port, settings = Some(ServerSettings(context.system)))
  override def receive: Receive = {
    case Http.Connected(remote, _) =>
      log.debug("Remote address {} connected", remote)
      sender ! Http.Register(context.actorOf(HttpConnectionHandler.props(remote, sender)))
  }
}

object HttpConnectionHandler {

  def props(remote: InetSocketAddress, connection: ActorRef): Props =
    Props(new HttpConnectionHandler(remote, connection))
}

class HttpConnectionHandler(remote: InetSocketAddress, connection: ActorRef) extends Actor with ActorLogging {

  // We need to know when the connection dies without sending a `Tcp.ConnectionClosed`
  context.watch(connection)

  def receive: Receive = {

    case HttpRequest(GET, uri, _, entity, _) => handleQuery(sender, GET, uri, entity)
    case HttpRequest(POST, uri, _, entity, _) => handleQuery(sender, POST, uri, entity)
    case _: Tcp.ConnectionClosed =>
      log.debug("Stopping, because connection for remote address {} closed", remote)
      context.stop(self)
    case Terminated(`connection`) =>
      log.debug("Stopping, because connection for remote address {} died", remote)
      context.stop(self)
  }

  private def handleQuery(sender: ActorRef, method : spray.http.HttpMethod, uri: spray.http.Uri, entity : spray.http.HttpEntity) = {
    val query = method match {
      case GET => uri.query.get("query")
      case POST => {
        entity.as[HttpForm] match {
          case Right(form) => {
            FormFieldExtractor(form).field("query").as[String] match {
              case Right(fieldval) => Option(fieldval)
              case Left(fieldval) => Option(None)
            }
          }
          case Left(form) => Option(None)
        }
      }
    }
    query match {
      case Some(query) => {
        import ExecutionContext.Implicits.global
        implicit val timeout = Timeout(2000, SECONDS)

        val future = context.actorSelection("/user/master") ? query recover {
          case _ => "Timeout error"
        }
        val result = Await.result(future, timeout.duration).asInstanceOf[String]
        log.debug("Sent to browser: " + result)
        sender ! HttpResponse(entity = result)
      }
      case None => sender ! HttpResponse(entity = "Empty query.")
    }
  }
}
