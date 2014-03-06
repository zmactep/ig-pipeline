package ru.biocad.ig.primer

import org.biojava3.core.sequence.DNASequence
import scala.util.{Failure, Try, Success}
import org.biojava3.core.sequence.transcription.TranscriptionEngine
import org.biojava3.core.sequence.compound.{DNACompoundSet, AminoAcidCompoundSet}
import org.ahocorasick.trie.Trie
import scala.collection.mutable


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
  def toDNA(strand: Option[String]): Option[String] = strand map {_.replace('U', 'T').replace('u', 't')}
  def getKmers(strand: Option[String], k: Int): Option[List[String]] = strand map{_.sliding(k, 1).toList}

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
}
