import akka.actor.{ActorPath, Props, ActorSystem}
import com.typesafe.config.ConfigFactory
import listeners.HttpListener
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

  val httpListener = system.actorOf(HttpListener.props("localhost", conf.getInt("ig-backend.http_port")), "http-listener")
  val m = system.actorOf(Props[Master], "master")

  val w1 = worker("master")
  val w2 = worker("master")
  val w3 = worker("master")

  system.awaitTermination()
  Console.println("Bye!")

  def worker(name: String) = system.actorOf(Props(
    new SimpleWorker(ActorPath.fromString(
      "akka://%s/user/%s".format(system.name, name)))))
}