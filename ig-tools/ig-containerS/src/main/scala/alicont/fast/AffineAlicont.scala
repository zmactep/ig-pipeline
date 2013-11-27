package alicont.fast

/**
 * Created with IntelliJ IDEA.
 * User: pavel
 * Date: 27.11.13
 * Time: 22:51
 */
abstract class AffineAlicont(maxheight : Int, query : String, gap_open : Int, gap_ext : Int, score_matrix : Array[Array[Int]])
  extends AbstractAlicont(maxheight, query, score_matrix) {

  protected val _gap_open = gap_open
  protected val _gap_ext  = gap_ext
  protected val _deletionMatrix  = new Matrix(query.size, maxheight)
  protected val _insertionMatrix = new Matrix(query.size, maxheight)

  def pop() : Unit = {
    val ls = _strings.pop()
    _scoreMatrix.move(-ls.size)
    _deletionMatrix.move(-ls.size)
    _insertionMatrix.move(-ls.size)
  }
}
