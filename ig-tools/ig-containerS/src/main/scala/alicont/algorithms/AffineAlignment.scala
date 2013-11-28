package alicont.algorithms

import alicont.Matrix

/**
 * Created with IntelliJ IDEA.
 * User: pavel
 * Date: 27.11.13
 * Time: 15:33
 */
trait AffineAlignment {
  def extendMatrix(s : String, query : String, gapOpen : Int, gapExtend : Int, score_matrix : Array[Array[Int]],
                   insertion_matrix : Matrix, deletion_matrix : Matrix, matrix : Matrix) : Unit
  def traceback(s : String, query : String, score_matrix : Array[Array[Int]],
                insertion_matrix : Matrix, deletion_matrix : Matrix, matrix : Matrix) : (Int, (String, String))
}
