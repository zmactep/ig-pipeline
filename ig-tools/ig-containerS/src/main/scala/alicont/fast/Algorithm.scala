package alicont.fast

/**
 * Created with IntelliJ IDEA.
 * User: mactep
 * Date: 11.11.13
 * Time: 16:19
 */
object Algorithm {
  trait AlgorithmStandard{
    def extendAlign(s : String, query : String,
                    gap : Int, score_matrix : Array[Array[Int]],
                    matrix : Matrix) : Unit
    def traceback(s : String, query : String,
                  gap : Int, score_matrix : Array[Array[Int]],
                  alignment : Matrix) : (Int, (String, String))
  }
  trait AlgorithmAffinyGap{
    def extendAlign(s : String, query : String,
                    gapOpen : Int, gapExtend : Int,
                    score_matrix : Array[Array[Int]], insertion_matrix : Array[Array[Int]],
                    deletion_matrix : Array[Array[Int]],
                    matrix : Matrix) : Unit
    def traceback(s : String, query : String,
                  gapOpen : Int, gapExtend : Int, score_matrix : Array[Array[Int]],
                  alignment : Matrix) : (Int, (String, String))
  }

  abstract class NeedlemanWunschBase extends AlgorithmStandard {
    def extendAlign(s : String, query : String,
        gap : Int, score_matrix : Array[Array[Int]],
        matrix : Matrix) : Unit
     = {
      if (matrix.height == 0) {
        initMatrix(matrix, query.size, gap)
      }
      (1 to s.size).foreach(i => {
        matrix.move(1)
        matrix.last(0) = matrix.pred(0) + gap
        (1 to query.size).foreach(j => {
          val score = score_matrix(s(i - 1))(query(j - 1))
          matrix.last(j) = (matrix.pred(j - 1) + score :: matrix.pred(j) + gap :: matrix.last(j - 1) + gap :: Nil).max
        })
      })
    }

    def traceback(s : String, query : String,
                  gap : Int, score_matrix : Array[Array[Int]],
                  alignment : Matrix) : (Int, (String, String)) = {
      var (score, i, j, result_s, result_q) = initTraceback(alignment, s, query)

      while (i != 0 || j != 0) {
        val cs : Char = if (i > 0) s(i - 1) else 0
        val cq : Char = if (j > 0) query(j - 1) else 0
        if (j == 0 || (i != 0 && alignment(i)(j) == alignment(i - 1)(j) + gap)) {
          i -= 1
          result_s.append(cs)
          result_q.append('-')
        } else if (i == 0 || (j != 0 && alignment(i)(j) == alignment(i)(j - 1) + gap)) {
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
      (score, (result_q.reverse.toString(), result_s.reverse.toString()))
    }
    protected def initMatrix(matrix : Matrix, n : Int, gap : Int) : Unit
    protected def initTraceback(alignment : Matrix, s : String, q : String)
    : (Int, Int, Int, StringBuilder, StringBuilder)
  }
  object NeedlemanWunsch extends NeedlemanWunschBase {
    protected def initMatrix(matrix : Matrix, n : Int, gap : Int) = {
      matrix.move(1)
      (0 to n).foreach(i => matrix.last(i) = i * gap)
    }

    protected def initTraceback(alignment : Matrix, s : String, q : String) = {
      (alignment.last.last, s.size, q.size, new StringBuilder(), new StringBuilder())
    }
  }

  object SmithWaterman extends AlgorithmStandard {
    def extendAlign(s : String, query : String,
                    gap : Int, score_matrix : Array[Array[Int]],
                    matrix : Matrix) : Unit
    = {
      if (matrix.height == 0) {
        matrix.move(1)
        (0 to query.size).foreach(i => matrix.last(i) = 0)
      }
      (1 to s.size).foreach(i => {
        matrix.move(1)
        matrix.last(0) = matrix.pred(0)
        (1 to query.size).foreach(j => {
          val score = score_matrix(s(i - 1))(query(j - 1))
          matrix.last(j) = (matrix.pred(j - 1) + score
            :: matrix.pred(j) + gap :: matrix.last(j - 1) + gap :: 0 :: Nil).max
        })
      })
    }

    def traceback(s : String, query : String,
                  gap : Int, score_matrix : Array[Array[Int]],
                  alignment : Matrix) : (Int, (String, String)) = {
      var (i, j) = (0, 0)
      var score = score_matrix(0)(0)
      val result_s = new StringBuilder()
      val result_q = new StringBuilder()

      for (k <- 0 to s.size; l <- 0 to query.size) {
        if(alignment(k)(l) > score) {
          score = alignment(k)(l)
          i = k
          j = l
        }
      }

      while (i != 0 && j != 0 && alignment(i)(j) != 0) {
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
        } else if (alignment(i)(j) != 0) {
          assert(false)
        }
      }

      (score, (result_q.reverse.toString(), result_s.reverse.toString()))
    }
  }
  object SemiGlobal extends NeedlemanWunschBase {

    protected def initMatrix(matrix : Matrix, n : Int, gap : Int) = {
      matrix.move(1)
      (0 to n).foreach(i => matrix.last(i) = 0)
    }

    protected def initTraceback(alignment : Matrix, s : String, query : String) = {
      var result_s = new StringBuilder()
      var result_q = new StringBuilder()
      var ind_s = 0
      var ind_q = 0

      var best_score_s_end = alignment(s.size)(0)
      var best_score_q_end = alignment(0)(query.size)

      for(i <- 1 to s.size) {
        if(alignment(i)(query.size) > best_score_q_end) {
          ind_s = i
          best_score_q_end = alignment(i)(query.size)
        }
      }

      for(j <- 1 to query.size) {
         if(alignment(s.size)(j) > best_score_s_end) {
           ind_q = j
           best_score_s_end = alignment(s.size)(j)
         }
      }

      var score = 0

      if(best_score_q_end > best_score_s_end) {
        score = best_score_q_end
        ind_q = query.size
        (ind_s to s.size - 1).foreach(i => {result_q.append('-'); result_s.append(s(i))})
      } else {
        score = best_score_s_end
        ind_s = s.size
        (ind_q to query.size - 1).foreach(i => {result_s.append('-'); result_q.append(query(i))})
      }

      (score, ind_s, ind_q, result_s.reverse, result_q.reverse)
    }
  }
}
