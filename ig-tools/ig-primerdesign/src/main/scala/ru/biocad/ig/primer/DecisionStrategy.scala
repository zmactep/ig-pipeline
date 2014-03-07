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