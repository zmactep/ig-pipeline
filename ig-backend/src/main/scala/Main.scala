import akka.actor.{ActorPath, Props, ActorSystem}
import com.typesafe.config.ConfigFactory
import java.net.InetSocketAddress
import java.sql.{Connection, DriverManager}
import listeners.{TcpListener, HttpListener}
import master.Master
import workers.SimpleWorker

/**
 * Created with IntelliJ IDEA.
 * User: Kos
 * Date: 25.09.13
 * Time: 8:32
 * To change this template use File | Settings | File Templates.
 */

object Main extends App {
  val system = ActorSystem("DB")
  val conf = ConfigFactory.load()

  val connectionString = "jdbc:mysql://" +
    conf.getString("ig-backend.db_host") + ":" +
    conf.getString("ig-backend.db_port") + "/" +
    conf.getString("ig-backend.db_name") + "?user=" +
    conf.getString("ig-backend.db_user") + "&password=" +
    conf.getString("ig-backend.db_password")

  var conn: Connection = null
  try {
    conn = DriverManager.getConnection(connectionString)
  } catch {
    case e: Exception => print("Error connecting to DB: " + e.getMessage); System.exit(1)
  }
  val httpListener = system.actorOf(HttpListener.props("localhost", conf.getInt("ig-backend.http_port")), "http-listener")
  val tcpListener = system.actorOf(TcpListener.props(new InetSocketAddress("localhost", conf.getInt("ig-backend.tcp_port"))), "tcp-listener")
  val m = system.actorOf(Master.props(conn), "master")

  val workers = for (i <- 0 until conf.getInt("ig-backend.workers")) yield worker("master")

  system.awaitTermination()
  Console.println("Bye!")

  def worker(name: String) = system.actorOf(Props(
    new SimpleWorker(ActorPath.fromString(
      "akka://%s/user/%s".format(system.name, name)))))
}