import alicont.slow.Alicont
import alicont.fast.{AlicontLocal, AlicontGlobal}
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
  }
}
