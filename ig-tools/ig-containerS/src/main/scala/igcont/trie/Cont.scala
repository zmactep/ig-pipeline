package igcont.trie

import scala.collection.mutable

/**
 * Created with IntelliJ IDEA.
 * User: mactep
 * Date: 04.10.13
 * Time: 12:50
 */
class Cont {
  private var _cont : mutable.ArrayBuffer[ContData] = mutable.ArrayBuffer.empty[ContData]

  def this(initnode : Node) = {
    this()
    _cont += new ContData(initnode)
  }

  def copyOf(cont : Cont) : Unit = {
    _cont = new mutable.ArrayBuffer[ContData](cont._cont.size)
    for(e <- cont._cont) {
      _cont += new ContData(e.node)
    }
  }

  def push(node : Node) : Unit = _cont += new ContData(node)

  def nodeOf(i : Int) : Node = _cont(i).node

  def dataOf(i : Int) : Any = _cont(i).data

  def setDataOf(i : Int, data : Any) = _cont(i).data = data

  def size : Int = _cont.size
}
