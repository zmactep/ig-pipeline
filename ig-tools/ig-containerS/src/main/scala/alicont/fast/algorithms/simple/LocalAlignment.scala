package alicont.fast.algorithms.simple

import alicont.fast.algorithms.SimpleAlignment
import alicont.fast.Matrix

/**
 * Created with IntelliJ IDEA.
 * User: pavel
 * Date: 27.11.13
 * Time: 15:40
 */
object LocalAlignment extends SimpleAlignment {
  def extendMatrix(s : String, query : String, gap : Int, score_matrix : Array[Array[Int]], matrix : Matrix) : Unit = {
    if (matrix.height == 0) {
      matrix.move(1)
      (0 to query.size).foreach(i => matrix.last(i) = 0)
    }
    (1 to s.size).foreach(i => {
      matrix.move(1)
      matrix.last(0) = matrix.pred(0)
      (1 to query.size).foreach(j => {
        val score = score_matrix(s(i - 1))(query(j - 1))
        matrix.last(j) = (matrix.pred(j - 1) + score :: matrix.pred(j) + gap :: matrix.last(j - 1) + gap :: 0 :: Nil).max
      })
    })
  }

  def traceback(s : String, query : String, gap : Int, score_matrix : Array[Array[Int]], matrix : Matrix) : (Int, (String, String)) = {
    var (i, j) = (s.size, query.size)
    val result_s = new StringBuilder()
    val result_q = new StringBuilder()

    var score = Int.MinValue
    for (it <- 0 to s.size; jt <- 0 to query.size) {
      if (score < matrix(it)(jt)) {
        score = matrix(it)(jt)
        i = it
        j = jt
      }
    }

    while (i != 0 && j != 0 && matrix(i)(j) != 0) {
      val cs : Char = if (i > 0) s(i - 1) else 0
      val cq : Char = if (j > 0) query(j - 1) else 0
      if (i != 0 && matrix(i)(j) == matrix(i - 1)(j) + gap) {
        i -= 1
        result_s.append(cs)
        result_q.append('-')
      } else if (j != 0 && matrix(i)(j) == matrix(i)(j - 1) + gap) {
        j -= 1
        result_s.append('-')
        result_q.append(cq)
      } else if (matrix(i)(j) == matrix(i - 1)(j - 1) + score_matrix(cs)(cq)) {
        i -= 1
        j -= 1
        result_s.append(cs)
        result_q.append(cq)
      } else if (matrix(i)(j) != 0) {
        assert(false)
      }
    }

    (score, (result_q.reverse.toString(), result_s.reverse.toString()))
  }
}
