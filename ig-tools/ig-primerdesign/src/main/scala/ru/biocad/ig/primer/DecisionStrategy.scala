package ru.biocad.ig.primer

/**
 * @author kfeodorov
 **/

sealed trait DecisionStrategy {
  /**
   * Choose next nucleotide. Sequence is prepended by dummyHead node with nucl="START". It points to the list of first nucleotides in Triquence
   * @param prev previous nucleotide. None if curr == dummyHead("START")
   * @param curr current nucleotide. "START" for dummyHead
   * @param next list of candidates
   * @param index current nucleotide index in protein string
   * @return next nucleotide from "next" @param
   */
  def next(prev: Option[String], curr: String, next: Seq[String], index: Int): Option[String]
}

/**
 * Simple strategy. Next nucleotide is chosen only according to it's weight
 * @param prob - map weight->nucleotide
 */
class SimpleProbabilityDecisionStrategy(prob: Map[String, Int]) extends DecisionStrategy{
  override def next(prev: Option[String], curr: String, next: Seq[String], index: Int): Option[String] = {
    val filteredProbs = prob.filter{e => next.contains(e._1)}
    if (filteredProbs.nonEmpty) Option(MiscUtils.sample(filteredProbs))
    else None
  }
}

/**
 * Strategy, based on codon frequency table.
 */
class CodonFrequencyDecisionStrategy extends DecisionStrategy{
  val equalsWeights = Map("A" -> 1000, "C" -> 1000, "G" -> 1000, "U" -> 1000)
  val total = 1581056
  val weights = Map( "UUU" -> 36738, "UCU" -> 13723, "UAU" -> 26077, "UGU" -> 8732,
                  "UUC" -> 26655, "UCC" -> 14004, "UAC" -> 19204, "UGC" -> 10920,
                  "UUA" -> 22000, "UCA" -> 12367, "UAA" -> 3217, "UGA" -> 1739,
                  "UUG" -> 22187, "UCG" -> 13729, "UAG" -> 426, "UGG" -> 23955,
                  "CUU" -> 18499, "CCU" -> 11563, "CAU" -> 21480, "CGU" -> 32051,
                  "CUC" -> 17321, "CCC" -> 9154, "CAC" -> 15490, "CGC" -> 33132,
                  "CUA" -> 6256, "CCA" -> 13508, "CAA" -> 23778, "CGA" -> 6132,
                  "CUG" -> 80431, "CCG" -> 34497, "CAG" -> 46720, "CGG" -> 9988,
                  "AUU" -> 47046, "ACU" -> 14447, "AAU" -> 29426, "AGU" -> 14953,
                  "AUC" -> 38317, "ACC" -> 35978, "AAC" -> 33762, "AGC" -> 25265,
                  "AUA" -> 8607, "ACA" -> 12981, "AAA" -> 52512, "AGA" -> 4660,
                  "AUG" -> 42762, "ACG" -> 23419, "AAG" -> 16855, "AGG" -> 3041,
                  "GUU" -> 29188, "GCU" -> 24683, "GAU" -> 50786, "GGU" -> 38610,
                  "GUC" -> 23896, "GCC" -> 39754, "GAC" -> 29331, "GGC" -> 44073,
                  "GUA" -> 17622, "GCA" -> 32557, "GAA" -> 60419, "GGA" -> 14209,
                  "GUG" -> 40384, "GCG" -> 50057, "GAG" -> 27971, "GGG" -> 17812)

  override def next(prev: Option[String], curr: String, next: Seq[String], index: Int): Option[String] = {
    if (index < 2) Option(MiscUtils.sample(equalsWeights.filterKeys(p => next.contains(p))))
    else {
       DnaUtils.toRNA(prev) match {
         case Some(pr) => val prefix = pr.toUpperCase + DnaUtils.toRNA(Some(curr)).get.toUpperCase
           val m = weights.filterKeys(_.startsWith(prefix)).filterKeys(p => next.contains(p.last.toString))
           if (m.nonEmpty) Option(MiscUtils.sample(m).last.toUpper.toString)
           else None
         case None => None
       }
     }
  }
}