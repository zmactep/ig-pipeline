package alicont

/**
 * Created with IntelliJ IDEA.
 * User: mactep
 * Date: 28.10.13
 * Time: 13:52
 */
object Algorithm {
  def needleman_wunsch(s : String, query : String,
                       gap : Int, score_matrix : Array[Array[Int]],
                       last_line : Array[Int] = null) : Array[Array[Int]] = {
    val matrix = Array.fill(s.size + 1, query.size + 1)(0)

    matrix(0) = if (last_line == null) (0 to query.size).map(_ * gap).toArray else last_line
    (1 to s.size).foreach(i => {
      matrix(i)(0) = matrix(i - 1)(0) + gap
      (1 to query.size).foreach(j => {
        val score = score_matrix(s(i - 1))(query(j - 1))
        matrix(i)(j) = (matrix(i - 1)(j - 1) + score :: matrix(i - 1)(j) + gap :: matrix(i)(j - 1) + gap :: Nil).max
      })
    })

    if (last_line == null) matrix else matrix.drop(1)
  }

  def traceback(s : String, query : String,
                gap : Int, score_matrix : Array[Array[Int]],
                alignment : Matrix) : (String, String) = {
    var (i, j) = (s.size, query.size)
    val result_s = new StringBuilder()
    val result_q = new StringBuilder()

    while (i != 0 || j != 0) {
      val cs : Char = if (i > 0) s(i - 1) else 0
      val cq : Char = if (j > 0) query(j - 1) else 0
      if (i != 0 && alignment(i)(j) == alignment(i - 1)(j) + gap) {
        i -= 1
        result_s.append(cs)
        result_q.append('-')
      } else if (j != 0 && alignment(i)(j) == alignment(i)(j - 1) + gap) {
        j -= 1
        result_s.append('-')
        result_q.append(cq)
      } else if (alignment(i)(j) == alignment(i - 1)(j - 1) + score_matrix(cs)(cq)) {
        i -= 1
        j -= 1
        result_s.append(cs)
        result_q.append(cq)
      } else {
        assert(false)
      }
    }

    (result_q.reverse.toString(), result_s.reverse.toString())
  }
}
