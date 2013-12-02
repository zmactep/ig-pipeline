package alicont.conts.affine

import alicont.algorithms.affine.GlobalAlignment
import alicont.conts.AffineAlicont

/**
 * Created with IntelliJ IDEA.
 * User: pavel
 * Date: 27.11.13
 * Time: 23:32
 */
class AlicontGlobal(maxheight : Int, query : String, gap_open : Double, gap_ext : Double, score_matrix : Array[Array[Double]])
  extends AffineAlicont(maxheight : Int, query : String, gap_open : Double, gap_ext : Double
    , score_matrix : Array[Array[Double]]) {

  def push(s : String) : Unit = {
    _strings.push(s)
    GlobalAlignment.extendMatrix(s, _query, _gap_open, _gap_ext, _score, _scoreMatrix, _insertionMatrix, _deletionMatrix)
  }

  def alignment() : (Double, (String, String)) = {
    GlobalAlignment.traceback(target, _query, _score, _scoreMatrix, _insertionMatrix, _deletionMatrix)
  }
}
