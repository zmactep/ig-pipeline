import org.scalatest.FunSpec
import ru.biocad.ig.primer._
import ru.biocad.ig.primer.ProteinTriequence
import scala.collection.mutable
import scala.util.Random

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
//        val ds = new SimpleProbabilityDecisionStrategy(Map("A" -> 5, "C" -> 5, "G" -> 5, "T" -> 5, "U" -> 5))
        val ds = new CodonFrequencyDecisionStrategy
        val hairpinMinLen = 8
        val pieces = 5
        val bodySize = protein.length * 3 / pieces - 5
        val overlap = 20
        val iterations = 1000
        var bestNucl: Option[String] = None
        var bestSplit: Option[List[Int]] = None
        var bestPrimersGcVariance: Double = 1000

        for ((nucl, i) <- triq.sampleStream(ds).filter{s => !containsHairPin(s, hairpinMinLen)}.zipWithIndex.take(iterations)) {
          if (i % 100 == 0) println(s"$i sequences tested")
          val indicies = DnaUtils.findOverlaps(nucl, DnaUtils.splitWithEqualGC(nucl, pieces, bodySize), overlap, bodySize)
          val gcVar: Double = DnaUtils.calcVarGc(DnaUtils.splitStrandToPrimers(nucl, indicies, overlap)).getOrElse(1000)
          if (gcVar < bestPrimersGcVariance) {
            println(s"New min GC variance per primer: $gcVar")
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

    it ("should support genetic algorithm") {
      import ru.biocad.ig.primer.ConversionUtils._
      import ru.biocad.ig.primer.MetricUtils._
      for (protein <- List("DIQMTQSPSSLSASVGDRVTITCKASQSVSSDVGWYQQKPGKAPKLLIYSGSNRYSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQDYYSPWTFGQGTKVEIK","EVQLVQSGAEVKKPGSSVKVSCKASGYTFTNYVINWVRQAPGQGLEWIGYNDPYNDVSKYNEKFKGRATITSDKSTSTAYMELSSLRSEDTAVYYCAKEGGGKYVYAMDSWGQGTTVTVSS")) {
        println(s"Test protein is $protein")
        val triq = ProteinTriequence(protein)
        val ds = new CodonFrequencyDecisionStrategy
        val primerSize = 21
        val genSize = 10
        val m = 10 //top m sequences are taken from generation
        val hairpinThreshold = 1000000000
        val overlapThreshold = 1000000000
        val iterations = 5
        val coef = (1., 1000.)
        implicit val ord: Ordering[(Double, String)] = Ordering.by(_._1)
        var queue = mutable.PriorityQueue[(Double, String)]()

        def getGeneration: Stream[(Double, String)] = triq.sampleStream(ds).
          filter{s => !isHairpin(s, hairpinThreshold)}.
          map{identity}.map{_.get}. //filter Nones
          map{s => (primersOverlapScore(Option(s), primerSize, coef).get , s)}.
          filter{_._1 < overlapThreshold}.
          take(genSize)

        //split each string into primers, accumulates set of all primers for each position in string,
        //and chooses random primer in each set to construct new string.
        // ex List(AAAA, BBBB, primerSize = 2) => List(AABB, BBAA) or List(AAAA, AABB) or List(BBAA, BBBB) and so on
        def mutate(strs: Seq[String]): List[String] = {
          val mm = new mutable.HashMap[Int, mutable.Set[String]] with mutable.MultiMap[Int, String]
          strs.foreach{s: String => augmentString(s).sliding(primerSize, primerSize).toList.zipWithIndex.foreach{p: (String, Int) => mm.addBinding(p._2, p._1)}}
          (0 until strs.size).toList.map{_ =>
            (for {(k, v) <- mm.toList} yield Random.shuffle(v.toList).head).mkString("")}
        }

        for (i <- 0 until iterations) {
          println(s"iteration $i")
          queue ++= getGeneration
          queue = mutable.PriorityQueue[(Double, String)]() ++
            mutate(queue.takeRight(m).toList.map{_._2}/*takes m seq with lowest score*/).
              filter{s => !isHairpin(Option(s), hairpinThreshold)}.
              map{s => (primersOverlapScore(Option(s), primerSize, coef).get , s)}
        }
        println(s"Result: ${queue.last}")
      }
    }
  }
}