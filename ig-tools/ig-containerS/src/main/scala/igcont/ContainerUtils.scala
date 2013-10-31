package igcont

import scala.util.Random

/**
 * Created with IntelliJ IDEA.
 * User: mactep
 * Date: 31.10.13
 * Time: 13:54
 */
object ContainerUtils {
  def warmup() : Unit = {
    def warmup_fun() = {
      val cont = new Container("ACGT", 'N')
      val rand = new Random()
      for(i <- 0 until 1000) {
        cont.push((0 until 100).map(_ => "ACGT"(rand.nextInt(4))).foldRight("")((c, s) => s + c), i.toString)
      }
    }

    (0 to 100).foreach(_ => warmup_fun())
  }

  def print_stats(cont : Container) : Unit = {
    val rc = Runtime.getRuntime

    printf("Current memory usage: %dMB\n", (rc.totalMemory() - rc.freeMemory()) / 1024 / 1024)
    printf("Compression: %.2f\n", cont.nodes.toFloat / cont.fullsize)
    printf("Records count: %d\n", cont.size)
  }
}
