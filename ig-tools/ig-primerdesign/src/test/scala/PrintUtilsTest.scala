import org.scalatest.FunSpec
import ru.biocad.ig.primer.{SimpleProbabilityDecisionStrategy, ProteinTriequence, DnaUtils, PrintUtils}

/**
 * @author kfeodorov
 * @since 22.02.14
 */
class PrintUtilsTest extends FunSpec {
  describe("PrintUtils") {
    it ("should print Sequence-based solution") {
      val protein = "EIVLTQSPSSVTASAGETVTINCKSSQSVLYSSNNKNYLAWYQQRPGQSPRLLIYWASTRESGVPDRFSGSGSGTDFTLTISSFQPEDAAVYYCQQGYSGVTFGQGTKVEIKRTVAAPSVLF"
      val codons = DnaUtils.proteinToCodonSets(Option(protein))
      val seq = codons.map(_.flatMap{codon => DnaUtils.flatten(Option(codon)).getOrElse(List())})
      val indicies = DnaUtils.findOverlaps(seq, 15, 30)
      PrintUtils.printOverlaps(seq, indicies, 15)
    }

    it ("should print Triquence-based solution") {
      import ru.biocad.ig.primer.ConversionUtils._
      import ru.biocad.ig.primer.DnaUtils._
      var protein: Sequence = "QLQLQESGGGLVQPGGSLRLSCAASGFTFTRYAMSWVRQAPGKGLEWVSAINSGGGSTYYADSVKGRFTISRDNAKNTVYLQLNSLKTEDMADVLVCRGGGTLGGR"
      Map(0 -> Set("E"), 1 -> Set("V"), 4 -> Set("V"), 29 -> Set("S"), 74 -> Set("S"), 78 -> Set("L"), 82 -> Set("M"), 86 -> Set("R"),
      87 -> Set("A"), 90 -> Set("T"), 94 -> Set("Y"), 95 -> Set("Y")).foreach{(subs: (Int, Set[String])) => protein = protein.addAlternativeAt(subs._1, subs._2)}

      val triq = ProteinTriequence(protein)
      val ds = new SimpleProbabilityDecisionStrategy(Map("A" -> 5, "C" -> 2000, "G" -> 2000, "T" -> 5, "U" -> 5))
      val nucl: Option[Sequence] = triq.sample(ds).map{s => println(s"$s straight"); println(s"${reverseComplementRNA(Option(s)).get} revcomp"); String2Sequence(s)}
      val indicies = DnaUtils.findOverlaps(nucl, 20, 40)
      PrintUtils.printOverlaps(nucl, indicies, 20)
    }

    it ("should work with Pasha's dataset :)") {
      import ru.biocad.ig.primer.ConversionUtils._
      import ru.biocad.ig.primer.DnaUtils._
      for (protein <- List("DIQMTQSPSSLSASVGDRVTITCKASQSVSSDVGWYQQKPGKAPKLLIYSGSNRYSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQDYYSPWTFGQGTKVEIK","EVQLVQSGAEVKKPGSSVKVSCKASGYTFTNYVINWVRQAPGQGLEWIGYNDPYNDVSKYNEKFKGRATITSDKSTSTAYMELSSLRSEDTAVYYCAKEGGGKYVYAMDSWGQGTTVTVSS")) {
        val triq = ProteinTriequence(protein)
        val ds = new SimpleProbabilityDecisionStrategy(Map("A" -> 5, "C" -> 5, "G" -> 5, "T" -> 5, "U" -> 5))
        println(protein)
        val nucl: Option[Sequence] = triq.sample(ds).map{s => println(s"$s straight"); println(s"${reverseComplementRNA(Option(s)).get} revcomp"); String2Sequence(s)}
        val indicies = DnaUtils.findOverlaps(nucl, 20, 40)
        PrintUtils.printOverlaps(nucl, indicies, 20)
      }
    }

    it ("should supports new splitWithEqualGC() primers splitter") {
      import ru.biocad.ig.primer.ConversionUtils._
      import ru.biocad.ig.primer.DnaUtils._
      for (protein <- List("DIQMTQSPSSLSASVGDRVTITCKASQSVSSDVGWYQQKPGKAPKLLIYSGSNRYSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQDYYSPWTFGQGTKVEIK","EVQLVQSGAEVKKPGSSVKVSCKASGYTFTNYVINWVRQAPGQGLEWIGYNDPYNDVSKYNEKFKGRATITSDKSTSTAYMELSSLRSEDTAVYYCAKEGGGKYVYAMDSWGQGTTVTVSS")) {
        println(s"Test protein is $protein")
        val triq = ProteinTriequence(protein)
        val ds = new SimpleProbabilityDecisionStrategy(Map("A" -> 5, "C" -> 5, "G" -> 5, "T" -> 5, "U" -> 5))
        val hairpinMinLen = 8
        val pieces = 10
        val bodySize = protein.length * 3 / pieces - 5
        val overlap = 10
        val iterations = 1000
        var bestNucl: Option[String] = None
        var bestSplit: Option[List[Int]] = None
        var bestPrimersGcVariance: Double = 1000

        for ((nucl, i) <- triq.sampleStream(ds).filter{s => !containsHairPin(s, hairpinMinLen)}.zipWithIndex.take(iterations)) {
          if (i % 100 == 0) println(s"$i sequences tested")
          val indicies = DnaUtils.findOverlaps(nucl, DnaUtils.splitWithEqualGC(nucl, pieces, bodySize), overlap, bodySize)
          val gcVar: Double = DnaUtils.calcVarGc(DnaUtils.splitStrandToPrimers(nucl, indicies, overlap)).getOrElse(1000)
          if (gcVar < bestPrimersGcVariance) {
            println(s"New best GC: $gcVar")
            bestPrimersGcVariance = gcVar
            bestNucl = nucl
            bestSplit = indicies
          }
        }
        println(s"Result:")
        println(s"${bestNucl.get} straight\n${reverseComplementRNA(bestNucl).get} revcomp")
        PrintUtils.printOverlaps(bestNucl map String2Sequence, bestSplit, overlap)
      }
    }
  }
}