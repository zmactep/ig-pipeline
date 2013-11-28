import alicont.fast.conts.simple.{AlicontSemiglobal, AlicontLocal, AlicontGlobal}
import alicont.Scoring
import igcont.anno.Anno
import org.scalatest.FlatSpec
import org.scalatest.matchers.ShouldMatchers

/**
 * Created with IntelliJ IDEA.
 * User: mactep
 * Date: 30.10.13
 * Time: 20:23
 */
class AlgoTest extends FlatSpec with ShouldMatchers {

  "Annotations" should "store right" in {
    val anno = new Anno(Array("Regions", "Genes", "Sites"))

    val rec = anno.createRecord("First", 20)

    val i = rec.setAnnotation(0, "Regions", "FR1")
    rec.setAnnotation(1, i._1, i._2)
    rec.setAnnotation(2, "Regions", "CDR1")
    rec.setAnnotation(3, "Genes", "V")

    rec.annotationOf(0)("Regions") should be ("FR1")
    rec.annotationOf(1)("Regions") should be ("FR1")
    rec.annotationOf(2)("Regions") should be ("CDR1")
    rec.annotationOf(3)("Genes") should be ("V")
  }

  "Alignment" should "count right score" in {
    val path : String = "../../data/BLOSUM62.txt"
    val a = new AlicontGlobal(11, "MEANLY", -5, Scoring.loadMatrix(path))
    a.push("PLE")
    a.push("ASANT")
    a.push("LY")

    val (score, _) = a.alignment()

    score should be (8)

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

    val e = new AlicontSemiglobal(11, "EASPTMEALYLY", -5, Scoring.loadMatrix(path))
    e.push("EAS")
    e.push("LY")
    val (score4, _) = e.alignment()
    score4 should be (15)
  }

  "Semiglobal" should "count right score" in {
    val path : String = "../../data/NUC_simple.txt"
    val a = new AlicontSemiglobal(11, "CAGCACTTGGATTCTCGG", -1, Scoring.loadMatrix(path))
    a.push("CAGCGTGG")
    val (score, _) = a.alignment()
    score should be (4)

    val b = new AlicontSemiglobal(11, "CAGCGAACACTTGGATTCTCGG", -1, Scoring.loadMatrix(path))
    b.push("CAGCGTGG")
    val (score2, _) = b.alignment()
    score2 should be (4)

    val c = new AlicontSemiglobal(11, "ACGTCAT", -1, Scoring.loadMatrix(path))
    c.push("TCATGCA")
    val (score3, _) = c.alignment()
    score3 should be (4)

    val d = new AlicontSemiglobal(11, "ACAGATA", -1, Scoring.loadMatrix(path))
    d.push("AGT")
    val (score4, _) = d.alignment()
    score4 should be (2)


    // FIXME: Killer sample! Semiglobal and local both has error when tring to align this.
    val e = new AlicontSemiglobal(40, "AAAAAAGAAAAAAAATGCCAAAAAAATTGG", -1, Scoring.loadMatrix(path))
    e.push("AAAAAAAAAAAAAAAAAAAAAATCTGTCGTGTTGGTTT")
    val (score5, _) = e.alignment()
  }
}
