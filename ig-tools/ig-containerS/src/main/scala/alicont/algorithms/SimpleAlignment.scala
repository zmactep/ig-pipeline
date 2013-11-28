package alicont.algorithms

import alicont.Matrix

/**
 * Created with IntelliJ IDEA.
 * User: pavel
 * Date: 27.11.13
 * Time: 15:33
 */
trait SimpleAlignment {
  def extendMatrix(s : String, query : String, gap : Int, score_matrix : Array[Array[Int]], matrix : Matrix) : Unit
  def traceback(s : String, query : String, gap : Int, score_matrix : Array[Array[Int]], matrix : Matrix) : (Int, (String, String))
}
