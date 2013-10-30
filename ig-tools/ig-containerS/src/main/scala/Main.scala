import igcont.Container

import scala.util.Random

/**
 * Created with IntelliJ IDEA.
 * User: mactep
 * Date: 04.10.13
 * Time: 13:10
 */

object Main{
  def load_test(i : Int) = {
    val cont = new Container("ACGT", 'N')
    val rand = new Random()

    val rc = Runtime.getRuntime
    val start = System.currentTimeMillis()

    for(i <- 0 until i*1000) {
      cont.push((0 until 400).map(_ => "ACGT"(rand.nextInt(4))).foldRight("")((c, s) => s + c), i.toString)
    }
    printf("%d) %.2f (%dMB)\n", i, (System.currentTimeMillis() - start) / 1000.0,
      (rc.totalMemory() - rc.freeMemory()) / 1024 / 1024)
  }

  def warmup() = {
    def warmup_fun() = {
      val cont = new Container("ACGT", 'N')
      val rand = new Random()
      for(i <- 0 until 1000) {
        cont.push((0 until 500).map(_ => "ACGT"(rand.nextInt(4))).foldRight("")((c, s) => s + c), i.toString)
      }
    }

    (0 to 100).foreach(_ => warmup_fun())
  }

  def main(args : Array[String]) = {
    warmup()
    (1 to 7).foreach(i => load_test(i))
  }
}