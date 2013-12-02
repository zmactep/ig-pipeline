package alicont.algorithms

import alicont.common.Matrix

/**
 * Created with IntelliJ IDEA.
 * User: pavel
 * Date: 27.11.13
 * Time: 15:33
 */
trait AffineAlignment {
  def extendMatrix(s : String, query : String, gapOpen : Double, gapExtend : Double, score_matrix : Array[Array[Double]],
                   insertion_matrix : Matrix, deletion_matrix : Matrix, matrix : Matrix) : Unit
  def traceback(s : String, query : String, score_matrix : Array[Array[Double]],
                insertion_matrix : Matrix, deletion_matrix : Matrix, matrix : Matrix) : (Double, (String, String))
}
