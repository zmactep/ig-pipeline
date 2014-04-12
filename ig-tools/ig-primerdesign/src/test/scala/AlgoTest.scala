import org.scalatest.FunSpec
import org.scalatest.Matchers._
import ru.biocad.ig.primer._
import ru.biocad.ig.primer.ProteinTriequence
import scala.collection.mutable
import scala.util.{Success, Try, Random}

/**
 * @author kfeodorov
 * @since 29.03.14
 */
class AlgoTest extends FunSpec {
  describe("Algorithms") {
    it ("Sequence-based strand generator with not fully overlapped primers") {
      val protein = "EIVLTQSPSSVTASAGETVTINCKSSQSVLYSSNNKNYLAWYQQRPGQSPRLLIYWASTRESGVPDRFSGSGSGTDFTLTISSFQPEDAAVYYCQQGYSGVTFGQGTKVEIKRTVAAPSVLF"
      val codons = DnaUtils.proteinToCodonSets(Option(protein))
      val seq = codons.map(_.flatMap{codon => DnaUtils.flatten(Option(codon)).getOrElse(List())})
      val indicies = DnaUtils.findOverlaps(seq, 15, 30)
      PrintUtils.printOverlaps(seq, indicies, 15)
    }

    it ("Triquence-based strand generator with not fully overlapped primers") {
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

    it ("Triquence-based strand generator with splitWithEqualGC() primers splitter") {
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

    it ("Genetic algorithm with Triquence-based strand generator and similarityScore()") {
      import ru.biocad.ig.primer.ConversionUtils._
      import ru.biocad.ig.primer.MetricUtils._
      for (protein <- List("DIQMTQSPSSLSASVGDRVTITCKASQSVSSDVGWYQQKPGKAPKLLIYSGSNRYSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQDYYSPWTFGQGTKVEIK","EVQLVQSGAEVKKPGSSVKVSCKASGYTFTNYVINWVRQAPGQGLEWIGYNDPYNDVSKYNEKFKGRATITSDKSTSTAYMELSSLRSEDTAVYYCAKEGGGKYVYAMDSWGQGTTVTVSS")) {
        println(s"Test protein is $protein")
        val triq = ProteinTriequence(protein)
        val ds = new CodonFrequencyDecisionStrategy
        val primerSize = 50
        val genSize = 30
        val m = 10 //top m sequences are taken from generation
        val hairpinThreshold = 1000000000
        val overlapThreshold = 1000000000
        val iterations = 5
        val coef = (1., 1.)
        implicit val ord: Ordering[(Double, String)] = Ordering.fromLessThan((t1, t2) => t1._1 > t2._1)
        val queue = mutable.PriorityQueue[(Double, String)]()

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
            (for {(k, v) <- mm.toList.sortBy(_._1)} yield Random.shuffle(v.toList).head).mkString("")}
        }

        for (i <- 0 until iterations) {
          Try(queue.dequeue()) match {
            case Success(best) =>
              println(s"iteration $i. Best (lowest) score is: ${best._1}")
              queue.enqueue(best)
            case _ => println(s"iteration $i. Queue is empty")
          }
          getGeneration.foreach(p => queue.enqueue(p))
          val sortedByScore = queue.dequeueAll

          (sortedByScore.take(m) ++ //take best seq (with lowest score) and store them as-is for the next iteration
            mutate(sortedByScore.take(m).toList.map{_._2}/*takes m seq with lowest score and mutate them*/).
              filter{s => !isHairpin(Option(s), hairpinThreshold)}.
              map{s => (primersOverlapScore(Option(s), primerSize, coef).get , s)}) foreach(p => queue.enqueue(p))
        }
        val ans = queue.dequeue()
        println(s"Best Score: ${ans._1}")
        println(s"Result is translated to ${DnaUtils.translate(Option(ans._2)).getOrElse("ERROR")}")

        val s: String = ans._2
        val rc: String = augmentString(DnaUtils.reverseComplementRNA(Option(s)).get).reverse
        val rcCut: String = augmentString(augmentString(rc).drop(primerSize/2)).dropRight(primerSize/2)
        def mergeLastShortPrimer(p: List[String]) = if (p.last.length < primerSize / 2) p.dropRight(2) :+ (p.init.last + p.last) else p
        val senseStrandPrimers = mergeLastShortPrimer(augmentString(s).sliding(primerSize, primerSize).toList)
        val antiSenseStrandPrimers = mergeLastShortPrimer(augmentString(rcCut).sliding(primerSize, primerSize).toList)

        println(senseStrandPrimers.mkString("|"))
        print(List.fill(primerSize/2)(" ").mkString(""))
        println(antiSenseStrandPrimers.mkString("|"))
        println("Primers:")
        val primers = senseStrandPrimers ++ antiSenseStrandPrimers.map{s => augmentString(s).reverse}
        for (i <- 0 until senseStrandPrimers.size - 1) {
          println(s"p${2*i} ${primers(i).replace("U", "U")}")
          println(s"p${2*i + 1} ${primers(i + senseStrandPrimers.size).replace("U", "U")}")
        }
        println(s"p${2*(senseStrandPrimers.size - 1)} ${senseStrandPrimers.last.replace("U", "U")}")
      }
    }
  }
}