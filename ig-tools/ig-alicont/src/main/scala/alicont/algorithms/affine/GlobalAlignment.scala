package alicont.algorithms.affine

import alicont.algorithms.AffineAlignment
import alicont.common.Matrix

/**
 * Created with IntelliJ IDEA.
 * User: Sergey Knyazev (sergey.n.knyazev@gmail.com)
 * Date: 28.11.13
 * Time: 15:12
 */
object GlobalAlignment extends AffineAlignment {
  def extendMatrix(s : String, query : String, gapOpen : Double, gapExtend : Double, score_matrix : Array[Array[Double]],
                   insertion_matrix : Matrix, deletion_matrix : Matrix, matrix : Matrix) : Unit = {
    if (matrix.height == 0) {
      insertion_matrix.move(1)
      deletion_matrix.move(1)
      matrix.move(1)

      matrix.last(0) = 0
      deletion_matrix.last(0) = Double.NegativeInfinity
      (1 to query.size).foreach(i => matrix.last(i) = gapOpen + gapExtend * i)
      (1 to query.size).foreach(i => deletion_matrix.last(i) = gapOpen + gapExtend * i)
      (0 to query.size).foreach(i => insertion_matrix.last(i) = Double.NegativeInfinity)

    }
    (1 to s.size).foreach(i => {

      insertion_matrix.move(1)
      deletion_matrix.move(1)
      matrix.move(1)

      insertion_matrix.last(0) = insertion_matrix.pred(0) + gapExtend
      deletion_matrix.last(0) = Double.NegativeInfinity
      matrix.last(0) = matrix.pred(0) + gapExtend

      (1 to query.size).foreach(j => {
        val score = score_matrix(s(i - 1))(query(j - 1))

        insertion_matrix.last(j) = Math.max(insertion_matrix.pred(j) + gapExtend,
          matrix.pred(j) + (gapOpen + gapExtend))
        deletion_matrix.last(j) = Math.max(deletion_matrix.last(j-1) + gapExtend,
          matrix.last(j-1) + (gapOpen + gapExtend)
        )
        matrix.last(j) = (matrix.pred(j - 1) + score :: insertion_matrix.last(j) :: deletion_matrix.last(j) :: Nil).max
      })
    })
  }

  def traceback(s : String, query : String, score_matrix : Array[Array[Double]],
                deletion_matrix : Matrix, insertion_matrix : Matrix, matrix : Matrix) : (Double, (String, String)) = {
    var (i, j) = (s.size, query.size)
    val result_s = new StringBuilder()
    val result_q = new StringBuilder()

    while (i != 0 || j != 0) {
      val cs : Char = if (i > 0) s(i - 1) else 0
      val cq : Char = if (j > 0) query(j - 1) else 0
      if (j == 0) {
        i -= 1
        result_s.append(cs)
        result_q.append('-')
      } else if (i == 0) {
        j -= 1
        result_s.append('-')
        result_q.append(cq)
      } else if (matrix(i)(j) == deletion_matrix(i)(j)) {
        i -= 1
        result_s.append(cs)
        result_q.append('-')
      } else if (matrix(i)(j) == insertion_matrix(i)(j)) {
        j -= 1
        result_s.append('-')
        result_q.append(cq)
      } else if (matrix(i)(j) == matrix(i - 1)(j - 1) + score_matrix(cs)(cq)) {
        i -= 1
        j -= 1
        result_s.append(cs)
        result_q.append(cq)
      } else {
        assert(false)
      }
    }

    (matrix.last.last, (result_q.reverse.toString(), result_s.reverse.toString()))
  }
}
