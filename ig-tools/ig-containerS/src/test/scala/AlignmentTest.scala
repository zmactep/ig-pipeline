import alicont.conts.simple.{AlicontSemiglobal, AlicontLocal, AlicontGlobal}
import alicont.Scoring
import org.scalatest.FlatSpec
import org.scalatest.matchers.ShouldMatchers

/**
 * Created with IntelliJ IDEA.
 * User: pavel
 * Date: 28.11.13
 * Time: 20:39
 */
class AlignmentTest extends FlatSpec with ShouldMatchers {
  val path : String = "../../data/BLOSUM62.txt"

  "SimpleGlobal" should "count right score" in {
    val a = new AlicontGlobal(11, "MEANLY", -5, Scoring.loadMatrix(path))
    a.push("PLE")
    a.push("ASANT")
    a.push("LY")

    val (score, _) = a.alignment()

    score should be (8)
  }

  "SimpleLocal" should "count right score" in {
    val b = new AlicontLocal(11, "MEANLYLY", -5, Scoring.loadMatrix(path))
    b.push("PLE")
    b.push("ASANT")
    b.push("LY")

    val (score1, _) = b.alignment()
    score1 should be (16)

    val c = new AlicontLocal(11, "MEANLY", -5, Scoring.loadMatrix(path))
    c.push("PLE")
    c.push("ASANT")
    c.push("LY")

    val (score2, _) = c.alignment()
    score2 should be (16)

    val d = new AlicontLocal(11, "MEANLY", -5, Scoring.loadMatrix(path))
    d.push("ME")
    d.push("A")
    d.push("L")
    d.push("Y")

    val (score3, _) = d.alignment()
    score3 should be (20)
  }

  it should "make right alignment" in {
    val path : String = "../../data/NUC_simple.txt"

    val f = new AlicontLocal(40, "AAAAAAGAAAAAAAATGCCAAAAAAATTGG", -1, Scoring.loadMatrix(path))
    f.push("AAAAAAAAAAAAAAAAAAAAAATCTGTCGTGTTGGTTT")
    val (_, (q6, t6)) = f.alignment()
    q6.replaceAll("-", "") should be ("AAAAAAGAAAAAAAATGCCAAAAAAATTGG")
    t6.replaceAll("-", "") should be ("AAAAAAAAAAAAAAAAAAAAAATCTGTCGTGTTGGTTT")
  }

  "SimpleSemiglobal" should "count right score" in {
    val e = new AlicontSemiglobal(11, "EASPTMEALYLY", -5, Scoring.loadMatrix(path))
    e.push("EAS")
    e.push("LY")
    val (score4, _) = e.alignment()
    score4 should be (15)
  }

  it should "count right score with NUC matrix" in {
    val path : String = "../../data/NUC_simple.txt"

    val a = new AlicontSemiglobal(11, "CAGCACTTGGATTCTCGG", -1, Scoring.loadMatrix(path))
    a.push("CAGCGTGG")
    val (score, (q, t)) = a.alignment()
    score should be (4)
    q.replaceAll("-", "") should be ("CAGCACTTGGATTCTCGG")
    t.replaceAll("-", "") should be ("CAGCGTGG")

    val b = new AlicontSemiglobal(11, "CAGCGAACACTTGGATTCTCGG", -1, Scoring.loadMatrix(path))
    b.push("CAGCGTGG")
    val (score2, (q2, t2)) = b.alignment()
    score2 should be (4)
    q2.replaceAll("-", "") should be ("CAGCGAACACTTGGATTCTCGG")
    t2.replaceAll("-", "") should be ("CAGCGTGG")

    val c = new AlicontSemiglobal(11, "ACGTCAT", -1, Scoring.loadMatrix(path))
    c.push("TCATGCA")
    val (score3, (q3, t3)) = c.alignment()
    score3 should be (4)
    q3.replaceAll("-", "") should be ("ACGTCAT")
    t3.replaceAll("-", "") should be ("TCATGCA")

    val d = new AlicontSemiglobal(11, "ACAGATA", -1, Scoring.loadMatrix(path))
    d.push("AGT")
    val (score4, (q4, t4)) = d.alignment()
    score4 should be (2)
    q4.replaceAll("-", "") should be ("ACAGATA")
    t4.replaceAll("-", "") should be ("AGT")
  }

  it should "make right alignment" in {
    val path : String = "../../data/NUC_simple.txt"

    val e = new AlicontSemiglobal(40, "AAAAAAGAAAAAAAATGCCAAAAAAATTGG", -1, Scoring.loadMatrix(path))
    e.push("AAAAAAAAAAAAAAAAAAAAAATCTGTCGTGTTGGTTT")
    val (_, (q5, t5)) = e.alignment()
    q5.replaceAll("-", "") should be ("AAAAAAGAAAAAAAATGCCAAAAAAATTGG")
    t5.replaceAll("-", "") should be ("AAAAAAAAAAAAAAAAAAAAAATCTGTCGTGTTGGTTT")
  }
}
