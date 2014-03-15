package ru.biocad.ig.primer

import scala.collection.mutable.ArrayBuffer

/**
 * @author kfeodorov
 * @since 15.03.14
 */
sealed trait Node
/**
 * Node of Trie representing nucleotide
 * @param nucl nucleotide in this node
 * @param next next nucleotides (in i + 1 th position) (with possible alternatives
 *             / synonymous substitutions) in reverse transcribed protein sequence
 */
final case class DataNode(nucl: String, next: ArrayBuffer[Node]) extends Node
object TerminalNode extends Node

object NodeUtils {
  def findNodeByNucl(nodes: Seq[Node], nucl: String): Option[Node] = nodes.find{case n: DataNode => n.nucl == nucl case _: TerminalNode.type => false}
  def findNodeByNuclWithChild(nodes: Seq[Node], nucl: String, nextNucl: String): Option[Node] = nodes.find {
    case n: DataNode => n.nucl == nucl && n.next.map{case s: DataNode => s.nucl case _:TerminalNode.type =>}.contains(nextNucl)
    case _: TerminalNode.type => false
  }
  def getChildren(node: Node) = node match {
    case d: DataNode => d.next
    case _: TerminalNode.type => Seq()
  }
}