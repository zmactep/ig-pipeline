import org.scalatest.FunSpec
import ru.biocad.ig.primer.{DnaUtils, PrintUtils}

/**
 * Created by Kos on 22.02.14.
 */
class PrintUtilsTest extends FunSpec {
  describe("PrintUtils") {
    it ("should just print :)") {
      val protein = "EIVLTQSPSSVTASAGETVTINCKSSQSVLYSSNNKNYLAWYQQRPGQSPRLLIYWASTRESGVPDRFSGSGSGTDFTLTISSFQPEDAAVYYCQQGYSGVTFGQGTKVEIKRTVAAPSVLF"
      val codons = DnaUtils.proteinToCodonSets(Option(protein))
      val seq = codons.map(_.flatMap{codon => DnaUtils.flattenCodons(Option(codon)).getOrElse(List())})
      val indicies = DnaUtils.findOverlaps(seq, 15, 30)
      PrintUtils.printOverlaps(seq.getOrElse(List[Set[String]]()), indicies.getOrElse(List()), 15)
    }
  }
}