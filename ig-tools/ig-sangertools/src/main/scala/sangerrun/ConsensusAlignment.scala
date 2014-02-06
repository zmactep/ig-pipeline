package sangerrun

import scala.collection.mutable

/**
 * Modified semiglobal alignment. Works on sequences of bags of symbols rather than individual symbols.
 * For two sequences, packages each symbol into its own bag, and runs semiglobal alignment on bags as if on regular strings.
 * The resulting traces are merged into a new sequence of bags, now every bag containing two symbols (one possibly
 * being a gap). Subsequent input sequences are aligned against the accumulate sequence-of-bags from previous iteration.
 * The final result, the one sequence assembled from input Sanger reads, is constructed as a consensus of symbols in
 * every bag. Ties in consensus for a position are resolved by taking one of the equally most frequent symbols.
 * To prevent gaps from appearing at ends of the resulting sequence, gap symbols do not participate in
 * consensus if every bag to the left or right up to sequence start or end also contains a gap.
 *
 * Example:
 * [W][A][T][L][R]
 * [T][D][R][S][W]
 * [L][R][A][W][A]
 * ->
 * [W][A][T][L][R][-][-][-]
 * [-][-][T][D][R][S][W][-]
 * [-][-][-][L][R][A][W][A]
 * ->
 * [W--][A--][TT-][LDL][RRR][-SA][-WW][--A]
 * ->
 * WATLRSWA
 * (gaps at edges ignored in calculating consensus)
 */
object ConsensusAlignment {
  private val MATCH_SCORE = 1
  private val MISMATCH_SCORE = -1
  private val GAP_SCORE = -1

  def merge(sequences: String*) : String = {
    if (sequences.length == 0) throw new IllegalArgumentException("Cannot merge 0 sequences.")

    var acc = new BagSequence(sequences.head)
    for (seq <- sequences.tail if seq.length > 0) {
      acc = align(acc, new BagSequence(seq))
    }

    acc.consensus.replace("-", "")
  }

  private def align(seq1: BagSequence, seq2: BagSequence) : BagSequence = {
    val mat = matrix(seq1, seq2)
    val (trace1, trace2) = backtrack(seq1, seq2, mat)
    trace1 merge trace2
  }

  private def score(b1: Bag, b2: Bag) : Int = if (b1 intersects b2) MATCH_SCORE else MISMATCH_SCORE

  private def matrix(seq1: BagSequence, seq2: BagSequence) : Array[Array[Int]] = {

    val mat = Array.ofDim[Int](seq1.length + 1, seq2.length + 1)
    for (i <- 0 to seq1.length) mat(i)(0) = 0
    for (i <- 0 to seq2.length) mat(0)(i) = 0

    for (row <- 1 to seq1.length; col <- 1 to seq2.length) {
      mat(row)(col) = List(
        mat(row - 1)(col - 1) + score(seq1(row - 1), seq2(col - 1)),
        mat(row - 1)(col) + GAP_SCORE,
        mat(row)(col - 1) + GAP_SCORE
      ).max
    }

    mat
  }

  private def backtrack(seq1: BagSequence, seq2: BagSequence, mat: Array[Array[Int]]) : (BagSequence, BagSequence) = {
    val trace1 = new BagSequence("")
    val trace2 = new BagSequence("")

    def prependBags(b1: Bag, b2: Bag) : Unit = { trace1.prepend(b1); trace2.prepend(b2) }
    def prependSequences(s1: BagSequence, s2: BagSequence) : Unit =
      (s1 zip s2).toSeq.reverse.foreach { case (b1, b2) => prependBags(b1, b2) }

    val best_row = (1 to seq1.length).maxBy(row => mat(row)(seq2.length))
    val best_col = (1 to seq2.length).maxBy(col => mat(seq1.length)(col))

    var (row, col) = List( (best_row, seq2.length), (seq1.length, best_col) ).maxBy{ case (r, b) => mat(r)(b) }

    if (row == seq1.length)
      prependSequences(new BagSequence("-" * (seq2.length - col)), seq2.subsequence(col, seq2.length))
    else // col must be == seq2.length
      prependSequences(seq1.subsequence(row, seq1.length), new BagSequence("-" * (seq1.length - row)))

    while (row > 0 && col > 0) {
      val upperleft = mat(row - 1)(col - 1) + score(seq1(row - 1), seq2(col - 1))
      val upper = mat(row - 1)(col) + GAP_SCORE
      val left = mat(row)(col - 1) + GAP_SCORE
      mat(row)(col) match {
        case `upperleft` =>
          prependBags(seq1(row - 1), seq2(col - 1))
          row -= 1
          col -= 1
        case `upper` =>
          prependBags(seq1(row - 1), new Bag('-'))
          row -= 1
        case `left` =>
          prependBags(new Bag('-'), seq2(col - 1))
          col -= 1
        case _ => throw new RuntimeException("Bug in backtracking.")
      }
    }

    if (row == 0)
      prependSequences(new BagSequence("-" * col), seq2.subsequence(0, col))
    else // col must be 0
      prependSequences(seq1.subsequence(0, row), new BagSequence("-" * row))

    (trace1, trace2)
  }

  /**
   * Sequence of bags for the modified semiglobal alignment.
   * */
  private class BagSequence (cont: Seq[Bag]) extends Iterable[Bag] {
    def this(init: String) = this(mutable.ArrayBuffer[Bag](init.map(b => new Bag(b)) : _*))
    private val _cont = mutable.ArrayBuffer[Bag](cont : _*)

    override def iterator: Iterator[Bag] = _cont.iterator

    def apply(index: Int) : Bag = _cont(index)
    def length : Int = _cont.length
    def prepend(bag: Bag) : Unit = _cont.prepend(bag)
    def subsequence(start: Int, end: Int) = new BagSequence(_cont.slice(start, end))
    def merge(other: BagSequence) : BagSequence = {
      if (other.length != length) throw new IllegalArgumentException("Sequence differs in length.")
      _cont.zip(other._cont).foldRight(new BagSequence("")){ case ((bag1, bag2), acc) => acc.prepend(bag1 ++ bag2); acc }
    }
    def consensus : String = {
      val left = _cont.takeWhile(b => b.contains('-')).map(b => b.consensus(ignoreGaps = true))
      val right = _cont.reverse.takeWhile(b => b.contains('-')).reverse.map(b => b.consensus(ignoreGaps = true))
      val center = _cont.slice(left.length, _cont.length - right.length).map(b => b.consensus(ignoreGaps = false))
      (left ++ center ++ right).mkString("")
    }
  }
  
  /**
   * Element of BagSequence, a collection of characters.
   * @param init initial character to put in the bag
   */
  private class Bag(init: Character*) {
    private val _cont = mutable.ListBuffer[Character](init : _*)

    def add(c: Character) : Unit = _cont += c
    def contains(c: Character) : Boolean = _cont.contains(c)
    def intersects(other: Bag) : Boolean = _cont.exists(c => other.contains(c))
    def ++(other: Bag) : Bag = new Bag(_cont ++ other._cont : _*)
    def consensus(ignoreGaps: Boolean = false) : Character = {
      var groups = _cont.groupBy(identity)
      if (ignoreGaps) groups = groups.filter{ case (s, g) => s != '-' }
      val max_size = groups.maxBy{ case (_, duplicates) => duplicates.size }._2.size
      val most_common = groups.filter{ case (char, duplicates) => duplicates.size == max_size }
      val not_ag = most_common.filter(g => !"AG".contains(g._1)) // trying to avoid selecting A and G since these may result in a stop codon
      if (not_ag.size > 0) not_ag.head._1 else most_common.head._1
    }
  }
}