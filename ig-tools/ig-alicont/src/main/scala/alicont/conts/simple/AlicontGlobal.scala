package alicont.conts.simple

import alicont.algorithms.simple.GlobalAlignment
import alicont.conts.SimpleAlicont

/**
 * Created with IntelliJ IDEA.
 * User: pavel
 * Date: 27.11.13
 * Time: 23:32
 */
class AlicontGlobal(maxheight : Int, query : String, gap : Double, score_matrix : Array[Array[Double]])
  extends SimpleAlicont(maxheight : Int, query : String, gap : Double, score_matrix : Array[Array[Double]]) {

  def push(s : String) : Unit = {
    _strings.push(s)
    GlobalAlignment.extendMatrix(s, _query, _gap, _score, _scoreMatrix)
  }

  def alignment() : (Double, (String, String)) = {
    GlobalAlignment.traceback(target, _query, _gap, _score, _scoreMatrix)
  }
}
