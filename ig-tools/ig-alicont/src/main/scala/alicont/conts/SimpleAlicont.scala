package alicont.conts


/**
 * Created with IntelliJ IDEA.
 * User: pavel
 * Date: 27.11.13
 * Time: 22:50
 */
abstract class SimpleAlicont(maxheight : Int, query : String, gap : Double, score_matrix : Array[Array[Double]])
  extends AbstractAlicont(maxheight : Int, query : String, score_matrix : Array[Array[Double]]) {

  protected val _gap = gap

  def pop() : Unit = {
    val ls = _strings.pop()
    _scoreMatrix.move(-ls.size)
  }
}
