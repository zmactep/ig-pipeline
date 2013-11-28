package alicont.conts.affine

import alicont.algorithms.affine.SemiglobalAlignment
import alicont.AffineAlicont

/**
 * Created with IntelliJ IDEA.
 * User: pavel
 * Date: 27.11.13
 * Time: 23:37
 */
class AlicontSemiglobal(maxheight : Int, query : String, gap_open : Int, gap_ext : Int, score_matrix : Array[Array[Int]])
  extends AffineAlicont(maxheight : Int, query : String, gap_open : Int, gap_ext : Int, score_matrix : Array[Array[Int]]) {

  def push(s : String) : Unit = {
    _strings.push(s)
    SemiglobalAlignment.extendMatrix(s, _query, _gap_open, _gap_ext, _score, _scoreMatrix, _insertionMatrix, _deletionMatrix)
  }

  def alignment() : (Int, (String, String)) = {
    SemiglobalAlignment.traceback(target, _query, _score, _scoreMatrix, _insertionMatrix, _deletionMatrix)
  }
}