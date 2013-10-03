package listeners

import java.net.InetSocketAddress
import akka.actor._
import akka.io.{Tcp, IO}
import scala.concurrent.{Await, ExecutionContext}
import akka.util.{ByteString, Timeout}
import akka.pattern.ask
import scala.concurrent.duration._

/**
 * Created with IntelliJ IDEA.
 * User: Kos
 * Date: 29.09.13
 * Time: 22:41
 * To change this template use File | Settings | File Templates.
 */
object TcpListener {

  def props(endpoint: InetSocketAddress): Props = Props(new TcpListener(endpoint))
}

class TcpListener(endpoint: InetSocketAddress) extends Actor with ActorLogging {

  import context.system

  IO(Tcp) ! Tcp.Bind(self, endpoint)

  override def receive: Receive = {
    case Tcp.Connected(remote, _) =>
      log.debug("Remote address {} connected", remote)
      sender ! Tcp.Register(context.actorOf(TcpConnectionHandler.props(remote, sender)))
  }
}

object TcpConnectionHandler {

  def props(remote: InetSocketAddress, connection: ActorRef): Props =
    Props(new TcpConnectionHandler(remote, connection))
}

class TcpConnectionHandler(remote: InetSocketAddress, connection: ActorRef) extends Actor with ActorLogging {

  // We need to know when the connection dies without sending a `Tcp.ConnectionClosed`
  context.watch(connection)

  def receive: Receive = {
    case Tcp.Received(data) =>
      Option(data) match {
        case Some(query) => {
          import ExecutionContext.Implicits.global
          implicit val timeout = Timeout(2000, SECONDS)

          val future = context.actorSelection("/user/master") ? query.utf8String recover {
            case _ => "Timeout error"
          }
          val result = Await.result(future, timeout.duration).asInstanceOf[String]
          log.debug("Sent to tcp endpoint: " + result)
          sender ! Tcp.Write(ByteString(result))
        }
        case None => sender !  Tcp.Write(ByteString("Empty query."))
      }
    case _: Tcp.ConnectionClosed =>
      log.debug("Stopping, because connection for remote address {} closed", remote)
      context.stop(self)
    case Terminated(`connection`) =>
      log.debug("Stopping, because connection for remote address {} died", remote)
      context.stop(self)
  }
}
