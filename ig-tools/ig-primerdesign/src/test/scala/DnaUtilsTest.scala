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
        DnaUtils.reverseComplementDNA(Option("CGGC")) should be (Some("GCCG"))
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

    it ("Split strand into pieces with least GC-content variance") {
      val bigStrand = "GAAGUCCAACUUGUCCAGAGCGGCGCCGAGGUGAAAAAGCCGGGCAGCAGCGUCAAAGUGAGUUGCAAAGCUUCAGGCUACACGUUCACCAAUUACGUGAUAAACUGGGUGAGACAAGCACCUGGACAAGGAUUAGAAUGGAUAGGCUACAAUGACCCGU"
      DnaUtils.splitWithEqualGC(Option(bigStrand), pieces = 4, minLen = 40).get.head._2 should be (Set(List(0, 40, 80, 120, 160)))
      DnaUtils.splitWithEqualGC(Option(bigStrand), pieces = 4, minLen = 39).get.values.flatten.toSet should be (Set(List(0, 40, 80, 120, 160), List(0, 39, 81, 120, 160), List(0, 40, 80, 119, 160), List(0, 39, 82, 121, 160), List(0, 41, 81, 121, 160), List(0, 39, 80, 119, 160), List(0, 41, 80, 119, 160), List(0, 39, 79, 121, 160), List(0, 40, 79, 120, 160), List(0, 42, 81, 120, 160), List(0, 43, 82, 121, 160), List(0, 40, 79, 119, 160), List(0, 39, 78, 118, 160), List(0, 39, 78, 119, 160), List(0, 41, 82, 121, 160), List(0, 39, 78, 120, 160), List(0, 40, 81, 120, 160), List(0, 39, 80, 120, 160), List(0, 42, 81, 121, 160), List(0, 39, 81, 121, 160), List(0, 40, 82, 121, 160), List(0, 41, 80, 121, 160), List(0, 40, 79, 118, 160), List(0, 39, 79, 118, 160), List(0, 40, 79, 121, 160), List(0, 42, 82, 121, 160), List(0, 41, 81, 120, 160), List(0, 39, 79, 120, 160), List(0, 39, 80, 121, 160), List(0, 41, 80, 120, 160), List(0, 40, 80, 121, 160), List(0, 39, 78, 117, 160), List(0, 40, 81, 121, 160), List(0, 39, 79, 119, 160), List(0, 39, 78, 121, 160))      )
      DnaUtils.splitWithEqualGC(Option(bigStrand), pieces = 5, minLen = 40) should be (None)
      DnaUtils.splitWithEqualGC(Option(bigStrand), pieces = 4, minLen = 41) should be (None)
    }

    it ("split to sense and antisense primers set") {
      DnaUtils.splitStrandToOverlappingPrimers(Option("AACGGCTT"), 4) should be (Some(List("AACG", "GCTT"), List("GCCG")))
    }
  }
}