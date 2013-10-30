package igcont.trie

/**
 * Created with IntelliJ IDEA.
 * User: mactep
 * Date: 04.10.13
 * Time: 12:05
 */
class Impl {
  private val _root : Node = new Node()
  private var _size : Int  = 1

  def insert(node : Node, symbol : Char) : Node = {
    if (!node.contains(symbol)) {
      node.set(new Node(node, symbol, _size))
      _size += 1
    }
    node.get(symbol).get
  }

  def root : Node = _root

  def size : Int = _size
}
