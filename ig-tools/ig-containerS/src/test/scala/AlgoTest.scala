import alicont.conts.simple
import alicont.conts.affine
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

  "Simple global aligner" should "should pass the test" in {
    val pathBlosum : String = "../../data/BLOSUM62.txt"
    val a = new simple.AlicontGlobal(11, "MEANLY", -5, Scoring.loadMatrix(pathBlosum))
    a.push("PLE")
    a.push("ASANT")
    a.push("LY")

    val (score, (t,q)) = a.alignment()

    score should be (8)
    q.replaceAll("-", "") should be ("PLEASANTLY")
    t.replaceAll("-", "") should be ("MEANLY")
  }

  "Simple local aligner" should "pass the test" in {
    val pathBlosum : String = "../../data/BLOSUM62.txt"
    val pathSimple : String = "../../data/NUC_simple.txt"

    val b = new simple.AlicontLocal(11, "MEANLYLY", -5, Scoring.loadMatrix(pathBlosum))
    b.push("PLE")
    b.push("ASANT")
    b.push("LY")

    val (score1, (t1,q1)) = b.alignment()
    score1 should be (16)
    q1.replaceAll("-", "") should be ("PLEASANTLY")
    t1.replaceAll("-", "") should be ("MEANLYLY")

    val c = new simple.AlicontLocal(11, "MEANLY", -5, Scoring.loadMatrix(pathBlosum))
    c.push("PLE")
    c.push("ASANT")
    c.push("LY")

    val (score2, (t2,q2)) = c.alignment()
    score2 should be (16)
    q2.replaceAll("-", "") should be ("PLEASANTLY")
    t2.replaceAll("-", "") should be ("MEANLY")


    val d = new simple.AlicontLocal(11, "MEANLY", -5, Scoring.loadMatrix(pathBlosum))
    d.push("ME")
    d.push("A")
    d.push("L")
    d.push("Y")

    val (score3, (t3,q3)) = d.alignment()
    score3 should be (20)
    q3.replaceAll("-", "") should be ("MEALY")
    t3.replaceAll("-", "") should be ("MEANLY")

    val f = new simple.AlicontLocal(40, "AAAAAAGAAAAAAAATGCCAAAAAAATTGG", -1, Scoring.loadMatrix(pathSimple))
    f.push("AAAAAAAAAAAAAAAAAAAAAATCTGTCGTGTTGGTTT")
    val (_, (q6, t6)) = f.alignment()
    q6.replaceAll("-", "") should be ("AAAAAAGAAAAAAAATGCCAAAAAAATTGG")
    t6.replaceAll("-", "") should be ("AAAAAAAAAAAAAAAAAAAAAATCTGTCGTGTTGGTTT")

  }

  "Simple semiglobal aligner" should "pass the test" in {

    val pathBlosum : String = "../../data/BLOSUM62.txt"
    val pathSimple : String = "../../data/NUC_simple.txt"

    val a = new simple.AlicontSemiglobal(11, "CAGCACTTGGATTCTCGG", -1, Scoring.loadMatrix(pathSimple))
    a.push("CAGCGTGG")
    val (score, (q, t)) = a.alignment()
    score should be (4)
    q.replaceAll("-", "") should be ("CAGCACTTGGATTCTCGG")
    t.replaceAll("-", "") should be ("CAGCGTGG")

    val b = new simple.AlicontSemiglobal(11, "CAGCGAACACTTGGATTCTCGG", -1, Scoring.loadMatrix(pathSimple))
    b.push("CAGCGTGG")
    val (score2, (q2, t2)) = b.alignment()
    score2 should be (4)
    q2.replaceAll("-", "") should be ("CAGCGAACACTTGGATTCTCGG")
    t2.replaceAll("-", "") should be ("CAGCGTGG")

    val c = new simple.AlicontSemiglobal(11, "ACGTCAT", -1, Scoring.loadMatrix(pathSimple))
    c.push("TCATGCA")
    val (score3, (q3, t3)) = c.alignment()
    score3 should be (4)
    q3.replaceAll("-", "") should be ("ACGTCAT")
    t3.replaceAll("-", "") should be ("TCATGCA")

    val d = new simple.AlicontSemiglobal(11, "ACAGATA", -1, Scoring.loadMatrix(pathSimple))
    d.push("AGT")
    val (score4, (q4, t4)) = d.alignment()
    score4 should be (2)
    q4.replaceAll("-", "") should be ("ACAGATA")
    t4.replaceAll("-", "") should be ("AGT")

    val e = new simple.AlicontSemiglobal(40, "AAAAAAGAAAAAAAATGCCAAAAAAATTGG", -1, Scoring.loadMatrix(pathSimple))
    e.push("AAAAAAAAAAAAAAAAAAAAAATCTGTCGTGTTGGTTT")
    val (_, (q5, t5)) = e.alignment()
    q5.replaceAll("-", "") should be ("AAAAAAGAAAAAAAATGCCAAAAAAATTGG")
    t5.replaceAll("-", "") should be ("AAAAAAAAAAAAAAAAAAAAAATCTGTCGTGTTGGTTT")

    val f = new simple.AlicontSemiglobal(11, "EASPTMEALYLY", -5, Scoring.loadMatrix(pathBlosum))
    f.push("EAS")
    f.push("LY")
    val (score6, (t6,q6)) = f.alignment()
    score6 should be (15)
    q6.replaceAll("-", "") should be ("EASLY")
    t6.replaceAll("-", "") should be ("EASPTMEALYLY")
  }

  "Global aligner with affine gaps" should "pass the test" in {
    val pathBlosum : String = "../../data/BLOSUM62.txt"

    val a = new affine.AlicontGlobal(11, "MEANLY", 0, -5, Scoring.loadMatrix(pathBlosum))
    a.push("PLE")
    a.push("ASANT")
    a.push("LY")

    val (score, (t,q)) = a.alignment()

    score should be (8)
    q.replaceAll("-", "") should be ("PLEASANTLY")
    t.replaceAll("-", "") should be ("MEANLY")

  }

  "Local aligner with affine gaps" should "pass the test" in {
    val pathBlosum : String = "../../data/BLOSUM62.txt"
    val pathSimple : String = "../../data/NUC_simple.txt"

    val b = new affine.AlicontLocal(11, "MEANLYLY", 0, -5, Scoring.loadMatrix(pathBlosum))
    b.push("PLE")
    b.push("ASANT")
    b.push("LY")

    val (score1, (t1,q1)) = b.alignment()
    score1 should be (16)
    q1.replaceAll("-", "") should be ("PLEASANTLY")
    t1.replaceAll("-", "") should be ("MEANLYLY")

    val c = new affine.AlicontLocal(11, "MEANLY", 0, -5, Scoring.loadMatrix(pathBlosum))
    c.push("PLE")
    c.push("ASANT")
    c.push("LY")

    val (score2, (t2,q2)) = c.alignment()
    score2 should be (16)
    q2.replaceAll("-", "") should be ("PLEASANTLY")
    t2.replaceAll("-", "") should be ("MEANLY")

    val d = new affine.AlicontLocal(11, "MEANLY", 0, -5, Scoring.loadMatrix(pathBlosum))
    d.push("ME")
    d.push("A")
    d.push("L")
    d.push("Y")

    val (score3, (t3,q3)) = d.alignment()
    score3 should be (20)
    q3.replaceAll("-", "") should be ("MEALY")
    t3.replaceAll("-", "") should be ("MEANLY")

    val f = new affine.AlicontLocal(40, "AAAAAAGAAAAAAAATGCCAAAAAAATTGG", 0, -1, Scoring.loadMatrix(pathSimple))
    f.push("AAAAAAAAAAAAAAAAAAAAAATCTGTCGTGTTGGTTT")
    val (_, (q6, t6)) = f.alignment()
    q6.replaceAll("-", "") should be ("AAAAAAGAAAAAAAATGCCAAAAAAATTGG")
    t6.replaceAll("-", "") should be ("AAAAAAAAAAAAAAAAAAAAAATCTGTCGTGTTGGTTT")


  }

  "Semiglobal aligner with affine gaps" should "should pass the test" in {
    val pathBlosum : String = "../../data/BLOSUM62.txt"
    val pathSimple : String = "../../data/NUC_simple.txt"

    val a = new affine.AlicontSemiglobal(11, "CAGCACTTGGATTCTCGG", 0, -1, Scoring.loadMatrix(pathSimple))
    a.push("CAGCGTGG")
    val (score, (q, t)) = a.alignment()
    score should be (4)
    q.replaceAll("-", "") should be ("CAGCACTTGGATTCTCGG")
    t.replaceAll("-", "") should be ("CAGCGTGG")

    val b = new affine.AlicontSemiglobal(11, "CAGCGAACACTTGGATTCTCGG", 0, -1, Scoring.loadMatrix(pathSimple))
    b.push("CAGCGTGG")
    val (score2, (q2, t2)) = b.alignment()
    score2 should be (4)
    q2.replaceAll("-", "") should be ("CAGCGAACACTTGGATTCTCGG")
    t2.replaceAll("-", "") should be ("CAGCGTGG")

    val c = new affine.AlicontSemiglobal(11, "ACGTCAT", 0, -1, Scoring.loadMatrix(pathSimple))
    c.push("TCATGCA")
    val (score3, (q3, t3)) = c.alignment()
    score3 should be (4)
    q3.replaceAll("-", "") should be ("ACGTCAT")
    t3.replaceAll("-", "") should be ("TCATGCA")

    val d = new affine.AlicontSemiglobal(11, "ACAGATA", 0, -1, Scoring.loadMatrix(pathSimple))
    d.push("AGT")
    val (score4, (q4, t4)) = d.alignment()
    score4 should be (2)
    q4.replaceAll("-", "") should be ("ACAGATA")
    t4.replaceAll("-", "") should be ("AGT")


    val e = new affine.AlicontSemiglobal(40, "AAAAAAGAAAAAAAATGCCAAAAAAATTGG", 0, -1, Scoring.loadMatrix(pathSimple))
    e.push("AAAAAAAAAAAAAAAAAAAAAATCTGTCGTGTTGGTTT")
    val (_, (q5, t5)) = e.alignment()
    q5.replaceAll("-", "") should be ("AAAAAAGAAAAAAAATGCCAAAAAAATTGG")
    t5.replaceAll("-", "") should be ("AAAAAAAAAAAAAAAAAAAAAATCTGTCGTGTTGGTTT")

    val f = new affine.AlicontSemiglobal(11, "EASPTMEALYLY", 0, -5, Scoring.loadMatrix(pathBlosum))
    f.push("EAS")
    f.push("LY")
    val (score6, (t6,q6)) = f.alignment()
    score6 should be (15)
    q6.replaceAll("-", "") should be ("EASLY")
    t6.replaceAll("-", "") should be ("EASPTMEALYLY")

  }
}
