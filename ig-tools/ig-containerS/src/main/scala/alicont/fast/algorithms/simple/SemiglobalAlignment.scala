package alicont.fast.algorithms.simple

import alicont.fast.algorithms.SimpleAlignment
import alicont.fast.Matrix

/**
 * Created with IntelliJ IDEA.
 * User: pavel
 * Date: 27.11.13
 * Time: 15:40
 */
object SemiglobalAlignment extends SimpleAlignment {
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
        matrix.last(j) = (matrix.pred(j - 1) + score :: matrix.pred(j) + gap :: matrix.last(j - 1) + gap :: Nil).max
      })
    })
  }

  def traceback(s : String, query : String, gap : Int, score_matrix : Array[Array[Int]], matrix : Matrix) : (Int, (String, String)) = {
    var (score, i, j) = prepareIJ(s.size, query.size, matrix)
    val result_s = new StringBuilder()
    val result_q = new StringBuilder()

    while (i != 0 || j != 0) {
      val cs : Char = if (i > 0) s(i - 1) else 0
      val cq : Char = if (j > 0) query(j - 1) else 0
      if (i == 0) {
        j -= 1
        result_s.append('-')
        result_q.append(cq)
      } else if (j == 0) {
        i -= 1
        result_s.append(cs)
        result_q.append('-')
      } else if (i != 0 && matrix(i)(j) == matrix(i - 1)(j) + gap) {
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
      } else {
        assert(false)
      }
    }

    (score, (result_q.reverse.toString(), result_s.reverse.toString()))
  }

  def prepareIJ(s : Int, q: Int, matrix : Matrix) : (Int, Int, Int) = {
    var (maxlastrow, maxlastcol) = (Int.MinValue, Int.MinValue)
    var (maxi, maxj) = (0, 0)
    for (jt <- 0 to q) {
      if (maxlastrow < matrix.last(jt)) {
        maxlastrow = matrix.last(jt)
        maxj = jt
      }
    }
    for (it <- 0 to s) {
      if (maxlastcol < matrix(it).last) {
        maxlastcol = matrix(it).last
        maxi = it
      }
    }

    if (maxlastrow >= maxlastcol) (maxlastrow, s, maxj) else (maxlastcol, maxi, q)
  }
}
