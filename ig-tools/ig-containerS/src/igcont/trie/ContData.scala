package igcont.trie

/**
 * Created with IntelliJ IDEA.
 * User: mactep
 * Date: 04.10.13
 * Time: 12:46
 */
class ContData(n : Node) {
  private val _node : Node = n
  var data  : Any  = null

  def node : Node = _node
}
