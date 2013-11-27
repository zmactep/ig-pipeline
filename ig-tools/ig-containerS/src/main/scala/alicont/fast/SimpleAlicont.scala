package alicont.fast

/**
 * Created with IntelliJ IDEA.
 * User: pavel
 * Date: 27.11.13
 * Time: 22:50
 */
abstract class SimpleAlicont(maxheight : Int, query : String, gap : Int, score_matrix : Array[Array[Int]])
  extends AbstractAlicont(maxheight : Int, query : String, score_matrix : Array[Array[Int]]) {

  protected val _gap = gap

  def pop() : Unit = {
    val ls = _strings.pop()
    _scoreMatrix.move(-ls.size)
  }
}
