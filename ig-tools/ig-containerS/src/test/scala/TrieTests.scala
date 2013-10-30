import igcont.trie.Trie
import org.scalatest.FlatSpec
import org.scalatest.matchers.ShouldMatchers

/**
 * Created with IntelliJ IDEA.
 * User: mactep
 * Date: 30.10.13
 * Time: 15:22
 */
class TrieTests extends FlatSpec with ShouldMatchers {

  "A Trie" should "keep right order of elements" in {
    val t = new Trie()
    val sa = Array("ACGG", "GCGT", "ACGT", "GTC")

    sa.foreach(s => {
      var k = 0
      s.foreach(c => {
        k = t.insert(k, c)
      })
    })

    t.symbolOf(1) should be ('A')
    t.symbolOf(2) should be ('C')
    t.symbolOf(3) should be ('G')
    t.symbolOf(4) should be ('G')
    t.symbolOf(5) should be ('G')
    t.symbolOf(6) should be ('C')
    t.symbolOf(7) should be ('G')
    t.symbolOf(8) should be ('T')
    t.symbolOf(9) should be ('T')
    t.symbolOf(10) should be ('T')
    t.symbolOf(11) should be ('C')
  }

  it should "has right dfs-order" in {
    val t = new Trie()
    val sa = Array("ACGG", "GCGT", "ACGT", "GTC")

    sa.foreach(s => {
      var k = 0
      s.foreach(c => {
        k = t.insert(k, c)
      })
    })

    var list_counter = 0
    var i = t.nextOf(0)
    while (i != 0) {
      if (t.isLeaf(i)) list_counter += 1
      i = t.nextOf(i)
    }

    list_counter should be (4)
  }
}
