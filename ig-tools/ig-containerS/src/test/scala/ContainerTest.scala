import org.scalatest._
import org.scalatest.matchers.ShouldMatchers


import alicont.Scoring
import igcont.Container
import scala.collection.mutable

/**
 * Created with IntelliJ IDEA.
 * User: mactep
 * Date: 30.10.13
 * Time: 15:25
 */
class ContainerTest extends FlatSpec with ShouldMatchers {

  "Container" should "have both name and handle addresation" in {
    val cont = new Container("ACGT", 'N')

    cont.push("ACGTAGCTACGATGCGACGACGACGAGGATGTTGGTTT", "Seq1")
    cont.push("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", "SeqA")
    cont.push("AAAAAAAAAAAAAAAAAAAAAATCTGTCGTGTTGGTTT", "Seq2")

    cont.seq("Seq2") should be(cont.seq(2))
  }

  it should "has correct annotations" in {
    val cont = new Container("ACGT", 'N', Array("A1", "A2", "A3"), 7)

    cont.push("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGT", "Seq1")
    cont.push("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", "SeqA")

    val rec = cont.record("SeqA")
    rec.setAnnotation(2, "A1", "42")
    rec.setAnnotation(3, "A1", "42")
    rec.setAnnotation(5, "A1", "38")
    rec.setAnnotation(3, "A2", ":)")

    val d = cont.data("SeqA").toArray
    d(2)._2("A1") should be ("42")
    d(3)._2("A1") should be ("42")
    d(5)._2("A1") should be ("38")
    d(3)._2("A2") should be (":)")
  }

  it should "search for alignments right" in {
    val path : String = "../../data/NUC4.4.txt"
    val cont = new Container("ACGT", 'N', Array("A1", "A2", "A3"), 7)

    cont.push("ACGTAGCTACGATGCGACGACGACGAGGATGTTGGTTT", "Seq1")
    cont.push("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", "SeqA")
    cont.push("AAAAAAAAAAAAAAAAAAAAAATCTGTCGTGTTGGTTT", "Seq2")

    val set = new mutable.TreeSet[String]()
    cont.alignment("AAAAAAGAAAAAAAATGCCAAAAAAATTGG", -5, Scoring.loadMatrix(path), 2).foreach(align => set += align.name)

    set should be (mutable.TreeSet[String]("SeqA", "Seq2"))

    cont.alignment("AAAAAAGAAAAAAAATGCCAAAAAAATTGG", -5, Scoring.loadMatrix(path), 0.57).head.name should be ("Seq2")
  }

  it should "search with specials" in {
    val cont = new Container("ACGT", 'N', Array(), 3)

    cont.push("GCTGGT", "Seq1")
    cont.push("AAAAAA", "SeqA")
    cont.push("ACTTGT", "Seq2")

    val set = new mutable.TreeSet[String]()
    cont.find("NCTNGT").foreach(tpl => set += tpl._1)

    set should be (mutable.TreeSet[String]("Seq1", "Seq2"))
  }
}
