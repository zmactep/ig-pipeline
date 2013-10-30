import _root_.scala.util.Random
import igcont.trie.Trie

/**
 * Created with IntelliJ IDEA.
 * User: mactep
 * Date: 30.10.13
 * Time: 15:22
 */
object TrieTests {
  def trie_test() = {
    println("*** TRIE TEST ***")

    val t = new Trie()
    val rand = new Random()
    val rc = Runtime.getRuntime

    val start = System.currentTimeMillis()
    (1 until 10000).foreach(i => {
      var key = 0
      (1 to 100).foreach(j => {
        key = t.insert(key, (65+rand.nextInt(4)).toChar)
        t.setDataOf(key, rand.nextInt(100))
      })
    })

    println((System.currentTimeMillis() - start) / 1000.0)

    val m = (rc.totalMemory() - rc.freeMemory()) / 1024 / 1024
    println(m)

    println(t.size)
  }

  def trie_iter_test() = {
    println("*** TRIE ITER TEST ***")

    val t = new Trie()
    val sa = Array("ACGG", "GCGT", "ACGT", "GTC")

    sa.foreach(s => {
      var k = 0
      s.foreach(c => {
        k = t.insert(k, c)
        print(k,"")
      })
      println()
    })

    var i = t.nextOf(0)
    while (i != 0) {
      print(i, "")
      i = t.nextOf(i)
    }
    println()
  }
}
