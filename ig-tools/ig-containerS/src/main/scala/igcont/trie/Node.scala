package igcont.trie

import scala.collection.mutable

/**
 * Created with IntelliJ IDEA.
 * User: mactep
 * Date: 04.10.13
 * Time: 11:47
 */
class Node(p: Node, s: Char, i: Int) {
  private val _parent   : Node            = p
  private val _symbol   : Char            = s
  private val _id       : Int             = i
  private val _children : mutable.HashMap[Char, Node] = new mutable.HashMap[Char, Node]()

  def this() = this(null, 0, 0)

  def parent : Node = _parent

  def symbol : Char = _symbol

  def id : Int = _id

  def keys : scala.collection.Set[Char] = _children.keySet

  def get(s : Char) : Option[Node] = _children.get(s)

  def set(n : Node) : Unit = _children.put(n.symbol, n)

  def contains(s : Char) : Boolean = _children.contains(s)

  def size : Int = _children.size
}
