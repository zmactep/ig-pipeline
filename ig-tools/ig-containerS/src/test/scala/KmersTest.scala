import _root_.scala.util.Random
import igcont.kmer.Counter

/**
 * Created with IntelliJ IDEA.
 * User: mactep
 * Date: 30.10.13
 * Time: 15:23
 */
object KmersTest {
  def kmers_test() = {
    println("*** KMERS TEST ***")

    val k = new Counter("01", '-', 3)
    val r = Array(0,1,2, 5, 3, 7, 6, 4)

    k.add(1, "0001011100", r)

    println(k.check("101", use_special = true))

    println(k.get("00-"))
    println(k.get("1-1"))
  }

  def kmers_nucleo_test(n : Int) = {
    println("*** KMERS NUCLEO TEST ***")

    val k = new Counter("ACGT", 'N', 5)
    val len = n - 4
    val r = 1 to len
    val rand = new Random()
    var s = ""
    (1 to n).foreach(i => s += "ACGT"(rand.nextInt(4)))
    println(s)

    k.add(1, s, r)

    println(k.get("GGNGT"))
    println(k.get("GGGGT"))
  }
}
