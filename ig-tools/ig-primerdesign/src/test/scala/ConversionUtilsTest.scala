import org.scalatest.FunSpec
import org.scalatest.Matchers._
import Predef.{augmentString => _, _}
import ru.biocad.ig.primer.DnaUtils.Sequence

/**
 * Created by Kos on 03.03.14.
 */
class ConversionUtilsTest extends FunSpec {
  describe("Conversions") {
    //Please note "Predef.{augmentString => _, _}" above - it hides ambiguous conversion augmentString()
    import ru.biocad.ig.primer.ConversionUtils._
    it ("implicitly converts string to Sequence") {
      List(Set("needed to trigger conversion")) ++ "F" should be (List(Set("needed to trigger conversion"), Set("F")))
      List(Set("needed to trigger conversion")) ++ "ACGT" should be (List(Set("needed to trigger conversion"), Set("A"), Set("C"), Set("G"), Set("T")))
    }

    it ("supports addition of alternatives to Sequence") {
      val seq: Sequence = "ACGT"
      seq.addAlternativeAt(0, Set("X")) should be (List(Set("A", "X"), Set("C"), Set("G"), Set("T")))
      seq.replaceAt(0, Set("X")) should be (List(Set("X"), Set("C"), Set("G"), Set("T")))
    }
  }
}
