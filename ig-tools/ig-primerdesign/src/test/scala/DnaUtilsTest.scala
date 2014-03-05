import org.scalatest.Matchers._
import org.scalatest.FunSpec
import ru.biocad.ig.primer.DnaUtils

/**
 * Created by Kos on 22.02.14.
 */
class DnaUtilsTest extends FunSpec {
  describe("DnaUtlis") {
    it ("should return correct reverse complement string") {
        DnaUtils.reverseComplementDNA(Option("A")) should be (Some("T"))
        DnaUtils.reverseComplementDNA(Option("t")) should be (Some("a"))
        DnaUtils.reverseComplementDNA(Option("acGT")) should be (Some("ACgt"))
        DnaUtils.reverseComplementDNA(Option("ACcgGT")) should be (Some("ACcgGT"))
    }

    it ("should accept N and -") {
      DnaUtils.reverseComplementDNA(Option("n")) should be (Some("n"))
      DnaUtils.reverseComplementDNA(Option("N")) should be (Some("N"))
      DnaUtils.reverseComplementDNA(Option("N---")) should be (Some("---N"))
    }

    it ("should ignore incorrect letters") {
      DnaUtils.reverseComplementDNA(Option("?")) should be (None)
      DnaUtils.reverseComplementDNA(Option("TAAAm")) should be (None)
    }

    it ("should reverse transcript aminoacid sequences") {
      DnaUtils.proteinToCodonSets(Option("FN")) should be (Some(List(Set("UUU", "UUC"), Set("AAC", "AAU"))))
    }

    it ("should flatten codons") {
      DnaUtils.flatten(Option(Set("UUU", "UUC", "UUA", "UUG"))) should be (Some(List(Set("U"), Set("U"), Set("U", "A", "C", "G"))))
      DnaUtils.flatten(Option(Set("U", "C", "A"))) should be (Some(List(Set("U", "A", "C"))))
      DnaUtils.flatten(Option(Set("U", "C", "AA"))) should be (Some(List(Set("U", "A", "C"), Set("A"))))
    }

    it ("should ignores wrong input") {
      val protein = "Extremely bad sequence"
      val codons = DnaUtils.proteinToCodonSets(Option(protein))
      val seq = codons.map(_.flatMap{codon => DnaUtils.flatten(Option(codon)).getOrElse(List())})
      seq should be (None)
    }

    it ("detect hairpins") {
      //http://en.wikipedia.org/wiki/File:Stem-loop.svg
      DnaUtils.containsHairPin(Option("ACGUGCCACGAUUCAACGUGGCACAG"), 8) should be (true)
      DnaUtils.containsHairPin(Option("ACGUGCCACGAUUCAACGUGGCACAG"), 9) should be (false)
    }

    it ("find duplicates") {
      DnaUtils.equals(Option(List(Set("A"), Set("C")), List(Set("A"), Set("C")))) should be (true)
      DnaUtils.equals(Option(List(Set("A"), Set("C")), List(Set("A"), Set("C"), Set("C")))) should be (true)
      DnaUtils.equals(Option(List(Set("A"), Set("C")), List(Set("A"), Set("C", "G"), Set("C")))) should be (false)
      DnaUtils.equals(Option(List(Set("A"), Set("C")), List(Set("A"), Set("C"), Set("C", "G")))) should be (true)
    }

    it ("find overlaps") {
      DnaUtils.findOverlaps(Option(List(Set("A"), Set("C"), Set("G"))), 1, 1) should be (Some(List(1)))
    }
  }
}