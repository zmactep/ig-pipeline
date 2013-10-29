package igcont.trie

/**
 * Created with IntelliJ IDEA.
 * User: mactep
 * Date: 04.10.13
 * Time: 14:18
 */
class Trie {
  private var _trie  : Impl       = new Impl()
  private val _cont  : Cont       = new Cont(_trie.root)
  private var _cache : Array[Int] = Array.empty[Int]

  def copyOf(trie : Trie) = {
    _trie = trie._trie
    _cont.copyOf(trie._cont)
    _cache = trie._cache
  }

  def insert(i : Int, symbol : Char) : Int = {
    val node = _trie.insert(_cont.nodeOf(i), symbol)
    if (node.id == _cont.size) {
      _cont.push(node)
    }
    node.id
  }

  def symbolOf(i : Int) : Char = _cont.nodeOf(i).symbol

  def setDataOf(i : Int, data : Any) : Unit = _cont.setDataOf(i, data)

  def dataOf(i : Int) : Any = _cont.dataOf(i)

  def parentOf(i : Int) : Option[Int] = {
    _cont.nodeOf(i).parent match {
      case null        => None
      case node : Node => Some(node.id)
    }
  }

  def nextOf(i : Int, symbol : Char) : Option[Int] = {
    _cont.nodeOf(i).get(symbol) match {
      case Some(n) => Some(n.id)
      case None    => None
    }
  }

  def nextOf(i : Int) : Int = {
    cache()
    _cache(i)
  }

  def keysOf(i : Int) : scala.collection.Set[Char] = _cont.nodeOf(i).keys

  def isFork(i : Int) : Boolean = _cont.nodeOf(i).size > 1

  def isLeaf(i : Int) : Boolean = _cont.nodeOf(i).size == 0

  def size : Int = _trie.size

  def cache() : Iterable[Int] = {
    def dfs(current: Int, last : Int) : Int = {
      _cache(last) = current
      var newlast : Int = current
      for (i <- _cont.nodeOf(current).keys) {
        newlast = dfs(_cont.nodeOf(current).get(i).get.id, newlast)
      }
      newlast
    }

    if (_cache.size != size) {
      _cache = Array.fill(size)(0)
      dfs(0, 0)
    }

    _cache
  }
}
