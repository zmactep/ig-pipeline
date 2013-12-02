package alicont.algorithms

import alicont.common.Matrix

/**
 * Created with IntelliJ IDEA.
 * User: pavel
 * Date: 27.11.13
 * Time: 15:33
 */
trait SimpleAlignment {
  def extendMatrix(s : String, query : String, gap : Double, score_matrix : Array[Array[Double]], matrix : Matrix) : Unit
  def traceback(s : String, query : String, gap : Double, score_matrix : Array[Array[Double]], matrix : Matrix)
  : (Double, (String, String))
}
