/*!
 * VisualEditor ContentEditable MWSyntaxHighlightNode class.
 *
 * @copyright 2011-2015 VisualEditor Team and others; see AUTHORS.txt
 * @license The MIT License (MIT); see LICENSE.txt
 */

/**
 * ContentEditable MediaWiki syntax highlight node.
 *
 * @class
 * @abstract
 *
 * @constructor
 */
ve.ce.MWSyntaxHighlightNode = function VeCeMWSyntaxHighlightNode() {
};

/* Inheritance */

OO.initClass( ve.ce.MWSyntaxHighlightNode );

/* Static Properties */

ve.ce.MWSyntaxHighlightNode.static.name = 'mwSyntaxHighlight';

ve.ce.MWSyntaxHighlightNode.static.primaryCommandName = 'syntaxhighlight';

/* Methods */

/** */
ve.ce.MWSyntaxHighlightNode.prototype.generateContents = function () {
	if ( !this.getModel().isLanguageSupported() ) {
		return $.Deferred().reject().promise();
	}
	// Parent method
	return ve.ce.MWExtensionNode.prototype.generateContents.apply( this, arguments );
};

/** */
ve.ce.MWSyntaxHighlightNode.prototype.onSetup = function () {
	// Parent method
	ve.ce.MWExtensionNode.prototype.onSetup.call( this );

	// DOM changes
	this.$element.addClass( 've-ce-mwSyntaxHighlightNode' );
};

/** */
ve.ce.MWSyntaxHighlightNode.prototype.getBoundingRect = function () {
	// HACK: Because nodes can overflow due to the pre tag, just use the
	// first rect (of the wrapper div) for placing the context.
	return this.rects[ 0 ];
};

/* Concrete subclasses */

ve.ce.MWBlockSyntaxHighlightNode = function VeCeMWBlockSyntaxHighlightNode() {
	// Parent method
	ve.ce.MWBlockExtensionNode.super.apply( this, arguments );

	// Mixin method
	ve.ce.MWSyntaxHighlightNode.call( this );
};

OO.inheritClass( ve.ce.MWBlockSyntaxHighlightNode, ve.ce.MWBlockExtensionNode );

OO.mixinClass( ve.ce.MWBlockSyntaxHighlightNode, ve.ce.MWSyntaxHighlightNode );

ve.ce.MWBlockSyntaxHighlightNode.static.name = 'mwBlockSyntaxHighlight';

ve.ce.MWInlineSyntaxHighlightNode = function VeCeMWInlineSyntaxHighlightNode() {
	// Parent method
	ve.ce.MWInlineExtensionNode.super.apply( this, arguments );

	// Mixin method
	ve.ce.MWSyntaxHighlightNode.call( this );
};

OO.inheritClass( ve.ce.MWInlineSyntaxHighlightNode, ve.ce.MWInlineExtensionNode );

OO.mixinClass( ve.ce.MWInlineSyntaxHighlightNode, ve.ce.MWSyntaxHighlightNode );

ve.ce.MWInlineSyntaxHighlightNode.static.name = 'mwInlineSyntaxHighlight';

/* Registration */

ve.ce.nodeFactory.register( ve.ce.MWBlockSyntaxHighlightNode );
ve.ce.nodeFactory.register( ve.ce.MWInlineSyntaxHighlightNode );
