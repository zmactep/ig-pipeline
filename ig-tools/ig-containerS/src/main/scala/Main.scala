import igcont.{ContainerUtils, Container}

import scala.collection.mutable.ArrayBuffer
import scala.util.Random

/**
 * Created with IntelliJ IDEA.
 * User: mactep
 * Date: 04.10.13
 * Time: 13:10
 */

object Main{
  private val strings = ArrayBuffer.empty[String]

  def generate() = {
    val rand = new Random()
    for(i <- 0 until 10000) {
      strings += (0 until 400).map(_ => "ACGT"(rand.nextInt(4))).foldRight("")((c, s) => s + c)
    }
  }

  def load_test(cont : Container, i : Int) = {
    val rc = Runtime.getRuntime
    val start = System.currentTimeMillis()

    (i * 1000 until (i+1) * 1000).foreach(j => cont.push(strings(j), j.toString))

    printf("%d) %.2f (%dMB)\n", i + 1, (System.currentTimeMillis() - start) / 1000.0,
      (rc.totalMemory() - rc.freeMemory()) / 1024 / 1024)
  }


  def main(args : Array[String]) = {
    println("Warmup")
    ContainerUtils.warmup()

    println("Generate")
    generate()

    println("Test")
    val cont = new Container("ACGT", 'N')
    (0 until 9).foreach(i => load_test(cont, i))
  }
}