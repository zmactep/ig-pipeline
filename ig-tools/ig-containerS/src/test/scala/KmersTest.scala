import igcont.kmer.Counter
import org.scalatest.FlatSpec
import org.scalatest.matchers.ShouldMatchers

/**
 * Created with IntelliJ IDEA.
 * User: mactep
 * Date: 30.10.13
 * Time: 15:23
 */
class KmersTest extends FlatSpec with ShouldMatchers {

  "K-mers" should "process special symbols correct" in {
    val k = new Counter("01", '-', 3)
    val r = Array(0, 1, 2, 5, 3, 7, 6, 4)

    k.add(1, "0001011100", r)

    val result = k.get("1-1").get._1

    result.head should be(5)
    result.last should be(7)
  }
}
