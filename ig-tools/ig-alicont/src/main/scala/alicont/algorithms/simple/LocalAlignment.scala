package alicont.algorithms.simple

import alicont.algorithms.SimpleAlignment
import alicont.common.Matrix

/**
 * Created with IntelliJ IDEA.
 * User: pavel
 * Date: 27.11.13
 * Time: 15:40
 */
object LocalAlignment extends SimpleAlignment {
  def extendMatrix(s : String, query : String, gap : Double, score_matrix : Array[Array[Double]], matrix : Matrix)
  : Unit = {
    if (matrix.height == 0) {
      matrix.move(1)
      (0 to query.size).foreach(i => matrix.last(i) = 0)
    }
    (1 to s.size).foreach(i => {
      matrix.move(1)
      matrix.last(0) = matrix.pred(0)
      (1 to query.size).foreach(j => {
        val score = score_matrix(s(i - 1))(query(j - 1))
        matrix.last(j) = (matrix.pred(j - 1) + score :: matrix.pred(j) + gap :: matrix.last(j - 1) + gap
          :: .0 :: Nil).max
      })
    })
  }

  def traceback(s : String, query : String, gap : Double, score_matrix : Array[Array[Double]], matrix : Matrix)
  : (Double, (String, String)) = {
    var (i, j) = (s.size, query.size)

    var score = Double.MinValue
    for (it <- 0 to s.size; jt <- 0 to query.size) {
      if (score < matrix(it)(jt)) {
        score = matrix(it)(jt)
        i = it
        j = jt
      }
    }

    var it = i
    var jt = j

    val end_s = new StringBuilder()
    val end_q = new StringBuilder()

    while(it != s.size || jt != query.size) {
      if (it == s.size) {
        end_s.append('-')
      } else {
        end_s.append(s(it))
        it += 1
      }
      if (jt == query.size){
        end_q.append('-')
      } else {
        end_q.append(query(jt))
        jt += 1
      }
    }

    val result_s = end_s.reverse
    val result_q = end_q.reverse

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
    while (i != 0 || j != 0) {
      if (i == 0) {
        result_s.append('-')
      } else {
        result_s.append(s(i - 1))
        i -= 1
      }
      if (j == 0){
        result_q.append('-')
      } else {
        result_q.append(query(j - 1))
        j -= 1
      }
    }

    (score, (result_q.reverse.toString(), result_s.reverse.toString()))
  }
}
