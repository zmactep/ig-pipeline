package ru.biocad.ig.primer

import org.biojava3.core.sequence.{RNASequence, DNASequence}
import scala.util.{Random, Failure, Try, Success}
import org.biojava3.core.sequence.transcription.TranscriptionEngine
import org.biojava3.core.sequence.compound.{DNACompoundSet, AminoAcidCompoundSet}
import org.ahocorasick.trie.Trie
import scala.collection.mutable
import scala.collection.immutable.SortedMap
import scala.annotation.tailrec
import scala.collection.mutable.ArrayBuffer


/**
 * Helper object with common routines
 */
object DnaUtils {
   type Sequence = List[Set[String]]

   val aminoAcid2CodonSet = MiscUtils.list2multimap({
    import scala.collection.JavaConversions._
    for {
      x <- TranscriptionEngine.getDefault.getTable.getCodons(DNACompoundSet.getDNACompoundSet, AminoAcidCompoundSet.getAminoAcidCompoundSet)
    } yield x.getAminoAcid.getShortName -> x.getShortName
  })

  val codon2aminoAcid = {
    import scala.collection.JavaConversions._
    for {
      x <- TranscriptionEngine.getDefault.getTable.getCodons(DNACompoundSet.getDNACompoundSet, AminoAcidCompoundSet.getAminoAcidCompoundSet)
    } yield x.getShortName -> x.getAminoAcid.getShortName
  }.toMap

  def reverseComplementDNA(strand: Option[String]): Option[String] = strand map { s => Try(new DNASequence(s).getReverseComplement.getSequenceAsString).getOrElse("") } filter(_.nonEmpty)
  def reverseComplementRNA(strand: Option[String]): Option[String] = strand map { s => Try(new RNASequence(s).getReverseComplement.getSequenceAsString).getOrElse("") } filter(_.nonEmpty)
  def toDNA(strand: Option[String]): Option[String] = strand map {_.replace('U', 'T').replace('u', 't')}
  def toRNA(strand: Option[String]): Option[String] = strand map {_.replace('T', 'U').replace('t', 'u')}
  def getKmers(strand: Option[String], k: Int): Option[List[String]] = strand map{_.sliding(k, 1).toList}
  def translate(rna: Option[String]): Option[String] = rna.map{_.sliding(3, 3).toList.map{codon2aminoAcid}.mkString("")}
  def getGC(strand: Option[String]): Option[Double] = strand.map{ str => str.toUpperCase.count{p: Char => 'C' == p || 'G' == p}.toDouble / str.length}

  /**
   * converts aminoacid string to codons list
   * @param protein string
   * @return list of sets of codons for i-th position
   */
  def proteinToCodonSets(protein: Option[String]): Option[Sequence] = protein match {
    case Some(p) => Try(for {
          acid <- p
        } yield aminoAcid2CodonSet(acid.toString)) match {
        case Success(s) => Option(s.toList)
        case Failure(_) => None
      }
    case None => None
  }

  /**
   * Converts list of codons => set of distinct nucleotides in i-th position.
   * Example: Option(Set("UUU", "UUC", "UUA", "UUG")) => Some(List(Set("U"), Set("U"), Set("U", "A", "C", "G")))
   * Actually, Strings in Set can be of any length, so this function is not limited to codons
   */
  def flatten(codon: Option[Set[String]]): Option[Sequence] = codon map { triplets =>
    //longest string in set
    val maxSize = triplets.foldLeft(0)((acc: Int, s: String) => List(acc, s.size).max)
    for {
      position <- 0 until maxSize
    } yield triplets.groupBy(s => Try(s.charAt(position)).getOrElse('*')).keySet.filter(_ != '*').map{_.toString}
  } map {_.toList}

  /**
   * Simple hairpin detector
   * @param strand DNA/RNA string
   * @param minPairedLen min len of complementary region
   * @return contains hair pin?
   */
  def containsHairPin(strand: Option[String], minPairedLen: Int): Boolean = strand match {
    case Some(s) => if (minPairedLen < 1) return false
      val trie = new Trie()
      getKmers(reverseComplementDNA(toDNA(Option(s))), minPairedLen) map {_.foreach{kmer => trie.addKeyword(kmer)} }
      !trie.parseText(toDNA(Option(s)).getOrElse("")).isEmpty
    case None => false
  }

  /**
   * Compare two Lists of Sets. i-th position in list contains a set of all possible nucleotides for this position in sequence
   * Two seqs are the same if for all positions in seq1 we have only one choice and it is the same as in seq2
   * @param seqs two sequences
   * @return true if sequences are the same
   */
  def equals(seqs: Option[(Sequence, Sequence)]): Boolean = seqs match {
    case Some((s1, s2)) => s1.zip(s2) forall {case (set1, set2) => set1.size == 1 && set1 == set2}
    case None => false
  }

  /**
   * Find unique regions in sequence that can be sticky ends (overlaps) in primers
   * @param strand input sequence
   * @param overlapSize sticky end size
   * @param interOverlapDistance non-overlapped part of primer. Non-edge primer will have length equals interOverlapDistance + 2 * overlapSize
   * @return indicies of sticky ends start
   */
  def findOverlaps(strand: Option[Sequence], overlapSize: Int, interOverlapDistance: Int): Option[List[Int]] = strand map { s =>
    def isUniq(seq2check: Sequence, others: List[Sequence]): Boolean = {
      if (others.isEmpty) true
      else others.forall(seq => !equals(Option((seq, seq2check))))
    }

    val ans = mutable.Buffer[Int]()
    val period = interOverlapDistance + overlapSize
    var pos = interOverlapDistance

    while (pos + interOverlapDistance < s.size) {
      //if current overlap is unique - choose it
        if (isUniq(seq2check = s.slice(pos, pos + overlapSize).toList, others = ans.map(i => s.slice(i, i + overlapSize).toList).toList)) {
          ans += pos
          pos += period
        } else pos += 1
    }
    ans.toList
  }

  /**
   * Splits strand into pieces of all random sizes (> minLen) and calculates gc-content variance for each variant of split.
   * @param strand strand to split
   * @param pieces # of pieces
   * @param minLen min len of a piece
   * @return gc content variance -> list of borders.
   */
  def splitWithEqualGC(strand: Option[String], pieces: Int, minLen: Int): Option[Map[Double, Set[List[Int]]]] = {
    /**
     * GC content as fraction
     * @param nom # 'G's and 'C's
     * @param denom total len
     */
    case class GC(nom: Int, denom: Int) {
      def double = nom.toDouble / denom
    }
    /**
     * Randomly generates borders.
     * ex: Border for first segment is chosen between (0 + minLen) and (size - (pieces - 1) * minLen)
     * @param size rightmost border position / total len of all segments
     * @param samples num of segments
     * @return list of segments
     */
    def genRandomBorders(size: Int, samples: Int): Set[List[Int]] = {
       def generate: List[Int] = {
         val ans = ArrayBuffer[Int](0)
           for (i <- 1 until pieces){
              val from = ans.last + minLen
              val to = size - minLen * (pieces - i)
              ans += from + Random.nextInt(to - from)
           }
           ans += size
         ans.toList
       }
      (for {i <- 0 until samples} yield generate).toSet// {
    }

    strand match {
      case Some(s) =>
        if (minLen * pieces > s.length) return None
        val cumSum = s.scanLeft(0)((acc: Int, c: Char) => if (List('G', 'C').contains(c.toUpper)) acc + 1 else acc)
        /**
         * gets GC content of segment of strand for O(1)
         * @param start index inclusive
         * @param stop index exclusive
         * @return GC content of segment
         */
        def getGC(start: Int, stop: Int) = GC(cumSum(stop) - cumSum(start), stop - start)

        //generate all possible borders
        val allBorders: Set[List[Int]] = genRandomBorders(s.size, 1000)
        //calc variance of GC for borders set and store in multimap gc_variance -> list of borders
        Some(MiscUtils.list2multimap(allBorders.toSeq.map{(list: List[Int]) => (MiscUtils.variance(for {
          segment: List[Int] <- list.sliding(2, 1).toList
        } yield getGC(segment.head, segment.last).double), list)}))
      case None => None
    }
  }

  /**
   * given borders, tries to build list of primers. Each primer has sticky end size of "overlap",
   * and has "body" size (primer - sticky ends) of at least bodySize.
   * @param dna DNA
   * @param segments splitWithEqualGC() result
   * @param overlap sticky end size
   * @param bodySize min primer body size
   * @return list of primer start indicies
   */
  def findOverlaps(dna: Option[String], segments: Option[Map[Double, Set[List[Int]]]], overlap: Int, bodySize: Int): Option[List[Int]] = {
    if (2 * overlap > bodySize) return None
    dna match {
      case Some(strand) =>
        segments match {
          case Some(seg) =>
            def isOverlapUniq(start: List[Int]): Boolean = start.forall{i => MiscUtils.countSubstring(strand, strand.substring(i, i + overlap)) == 1}
            val sMap = SortedMap[Double, Set[List[Int]]]() ++ seg.filter(_._2.nonEmpty)
            for (entry <- sMap) {
              for (borders: List[Int] <- entry._2) {
                val shifted = borders.tail.map{_ - overlap / 2}.init
                if (isOverlapUniq(shifted)) return Some(shifted)
              }
            }
            None
          case None => None
        }
      case None => None
    }
  }

  def splitStrandToPrimers(strand: Option[String], indicies: Option[List[Int]], overlap: Int): Option[List[String]] =
    strand.flatMap { s: String => indicies.map {idx: List[Int] => (0 :: idx).zip{idx :+ s.size}}.map{_.map{range: (Int, Int) => s.substring(range._1, math.min(range._2 + overlap, s.size))}} }

  def calcVarGc(strands: Option[List[String]]): Option[Double] = strands.map{list => MiscUtils.variance(list.map{s => getGC(Option(s)).get})}

  /**
   * splits strand to overlapped primers
   * @param strand
   * @param size primer len
   * @return -------------- =>
   *         -- --- --- --- sense
   *          --- --- ---   antisense
   *          returns (sense, antisense)
   */
  def splitStrandToOverlappingPrimers(strand: Option[String], size: Int): Option[(List[String], List[String])] = strand.map {
    s =>
      val primers = s.sliding(size, size / 2).toList
      val sense: List[String] = primers.zipWithIndex.filter{_._2 % 2 == 0}.unzip._1
      val antiSense: List[String] = primers.zipWithIndex.filter{_._2 % 2 == 1}.unzip._1.map{str: String => reverseComplementDNA(Option(str)).get}
      (sense, antiSense)
  }
}
