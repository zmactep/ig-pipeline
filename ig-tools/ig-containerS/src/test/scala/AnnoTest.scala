import igcont.anno.Anno
import org.scalatest.FlatSpec
import org.scalatest.matchers.ShouldMatchers

/**
 * Created with IntelliJ IDEA.
 * User: mactep
 * Date: 30.10.13
 * Time: 20:23
 */
class AnnoTest extends FlatSpec with ShouldMatchers {

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
}
