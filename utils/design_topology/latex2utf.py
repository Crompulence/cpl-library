#!/usr/bin/env python
# -*- coding: utf8 -*-

# latex2utf.py
# Copyright (C) 2008 Alexander Rodin                                        
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from re import sub
from sys import stdin

class Command:
    def __init__(self, is_callable = False, replacement = '', args_count = 0):
        self.is_callable = is_callable
        self.replacement = replacement
        self.args_count = args_count

mod_symbols = {
        u'0': '⁰',
        u'1': '¹',
        u'2': '²',
        u'3': '³',
        u'4': '⁴',
        u'5': '⁵',
        u'6': '⁶',
        u'7': '⁷',
        u'8': '⁸',
        u'9': '⁹',
        u'+': '⁺',
        u'-': '⁻',
        u'=': '⁼',
        u'n': 'ⁿ',
        u'(': '⁽',
        u')': '⁾',
        u'○': '°',
}

sub_symbols = {
        u'0': '₀',
        u'1': '₁',
        u'2': '₂',
        u'3': '₃',
        u'4': '₄',
        u'5': '₅',
        u'6': '₆',
        u'7': '₇',
        u'8': '₈',
        u'9': '₉',
        u'+': '₊',
        u'-': '₋',
        u'=': '₌',
        u'(': '₍',
        u')': '₎',
        u'a': 'ₐ',
        u'e': 'ₑ'
}

def is_str_in(str, a):
    s = str.decode('utf-8')
    for c in s:
        if not c in a:
            return False
    return True

def cmd_mod(arg):
    if is_str_in(arg.replace(' ', ''), mod_symbols.keys()):
        res = ''
        for c in arg.replace(' ', '').decode('utf-8'):
            res += mod_symbols[c]
        return res
    else:
        return '^(' + arg + ')'

    return '^(' + arg + ')'

def cmd_sub(arg):
    if is_str_in(arg.replace(' ', ''), sub_symbols.keys()):
        res = ''
        for c in arg.replace(' ', '').decode('utf-8'):
            res += sub_symbols[c]
        return res
    else:
        return '_(' + arg + ')'

    return '_(' + arg + ')'

def cmd_frac(args):
    a = args[0]
    b = args[1]
    if a == '1' and b == '4':
        return '¼'
    if a == '1' and b == '2':
        return '½'
    if a == '3' and b == '4':
        return '¾'
    if a == '1' and b == '3':
        return '⅓'
    if a == '2' and b == '3':
        return '⅔'
    if a == '1' and b == '5':
        return '⅕'
    if a == '2' and b == '5':
        return '⅖'
    if a == '3' and b == '5':
        return '⅗'
    if a == '4' and b == '5':
        return '⅘'
    if a == '1' and b == '6':
        return '⅙'
    if a == '5' and b == '6':
        return '⅚'
    if a == '1' and b == '8':
        return '⅛'
    if a == '5' and b == '8':
        return '⅝'
    if a == '7' and b == '8':
        return '⅞'
    if is_str_in(a.replace(' ', ''), mod_symbols.keys()) and is_str_in(b.replace(' ', ''), sub_symbols.keys()):
        return cmd_mod(a) + '⁄' + cmd_sub(b)
    return '(' + a + ')⁄(' + b + ')'

commands = {
        'alpha':        Command(False, 'α'),
        'beta':         Command(False, 'β'),
        'Gamma':        Command(False, 'Γ'),
        'gamma':        Command(False, 'γ'),
        'Delta':        Command(False, 'Δ'),
        'delta':        Command(False, 'δ'),
        'epsilon':      Command(False, 'ϵ'),
        'varepsilon':   Command(False, 'ε'),
        'zeta':         Command(False, 'ζ'),
        'eta':          Command(False, 'η'),
        'Theta':        Command(False, 'Θ'),
        'theta':        Command(False, 'θ'),
        'vartheta':     Command(False, 'ϑ'),
        'iota':         Command(False, 'ι'),
        'kappa':        Command(False, 'κ'),
        'Lambda':       Command(False, 'Λ'),
        'lambda':       Command(False, 'λ'),
        'mu':           Command(False, 'μ'),
        'nu':           Command(False, 'ν'),
        'Xi':           Command(False, 'Ξ'),
        'xi':           Command(False, 'ξ'),
        'Pi':           Command(False, 'Π'),
        'pi':           Command(False, 'π'),
        'varpi':        Command(False, 'ϖ'),
        'rho':          Command(False, 'ρ'),
        'varrho':       Command(False, 'ϱ'),
        'Sigma':        Command(False, 'Σ'),
        'sigma':        Command(False, 'σ'),
        'varsigma':     Command(False, 'ς'),
        'tau':          Command(False, 'τ'),
        'Upsilon':      Command(False, 'Υ'),
        'upsilon':      Command(False, 'υ'),
        'Phi':          Command(False, 'Φ'),
        'phi':          Command(False, 'ϕ'),
        'varphi':       Command(False, 'φ'),
        'chi':          Command(False, 'χ'),
        'Psi':          Command(False, 'Ψ'),
        'psi':          Command(False, 'ψ'),
        'Omega':        Command(False, 'Ω'),
        'omega':        Command(False, 'ω'),
        'int':          Command(False, '∫'),
        'cdot':         Command(False, '⋅'),
        'pm':           Command(False, ' ± '),
        'mp':           Command(False, ' ∓ '),
        'lor':          Command(False, ' ∨ '),
        'land':         Command(False, ' ∧ '),
        'le':           Command(False, ' ≤ '),
        'ge':           Command(False, ' ≥ '),
        'equiv':        Command(False, ' ≡ '),
        'sim':          Command(False, ' ∼ '),
        'parallel':     Command(False, '∥ '),
        'perp':         Command(False, ' ⊥ '),
        'infty':        Command(False, '∞'),
        'times':        Command(False, '⨯'),
        'll':           Command(False, ' ≪ '),
        'gg':           Command(False, ' ≫ '),
        'simeq':        Command(False, ' ≃ '),
        'approx':       Command(False, ' ≈ '),
        'neq':          Command(False, ' ≠ '),
        'angle':        Command(False, '∠'),
        'triangle':     Command(False, '△'),
        'sum':          Command(False, '∑'),
        'Rightarrow':   Command(False, ' ⇒ '),
        'Leftrightarrow':Command(False, ' ⇔ '),
        'wedge':        Command(False, '∧'),
        'vee':          Command(False, '∨'),
        'neg':          Command(False, '¬'),
        'forall':       Command(False, '∀'),
        'exists':       Command(False, '∃'),
        'varnothing':   Command(False, '∅'),
        'in':           Command(False, ' ∈ '),
        'notin':        Command(False, ' ∉ '),
        'subseteq':     Command(False, ' ⊆ '),
        'subset':       Command(False, ' ⊂ '),
        'cup':          Command(False, ' ∪ '),
        'cap':          Command(False, ' ⋂ '),
        'to':           Command(False, ' → '),
        'mapsto':       Command(False, ' ↦ '),
        'prod':         Command(False, '∏'),
        'circ':         Command(False, '○'),
        'sin':          Command(False, ' sin '),
        'cos':          Command(False, ' cos '),
        'tan':          Command(False, ' tan '),
        'ctab':         Command(False, ' ctan '),
        'asin':         Command(False, ' asin '),
        'acos':         Command(False, ' acos '),
        'atan':         Command(False, ' atan '),
        'actan':        Command(False, ' actan '),
        'log':          Command(False, ' log '),
        'ln':           Command(False, ' ln '),
        'lg':           Command(False, ' lg '),

        'mathbb':       Command(True, lambda args: { 'N': 'ℕ', 'Z': 'ℤ', 'Q': 'ℚ', 'R': 'ℝ', 'C': 'ℂ' }[args[0]], 1),
        'vec':          Command(True, lambda args: args[0] + '⃗', 1),
        'sqrt':         Command(True, lambda args: '√(' + args[0] + ')', 1),
        'frac':         Command(True, cmd_frac, 2)
}

def latex2utf(src):
    res = ''

    i = 0
    while i < len(src):
        if src[i] == '\\' and len(src) > i + 1:
            if src[i + 1] == '\\':
                res += '\\'
                i += 2
                continue
            cmd = ''
            for c in src[i + 1:]:
                if c.isalpha():
                    cmd += c
                else:
                    break
            i += len(cmd)
            if commands.__contains__(cmd):
                f = commands[cmd]
                if f.is_callable:
                    args = []
                    if f.args_count > 0:
                        for j in range(0, f.args_count):
                            i += 1
                            if i >= len(src):
                                return ""
                            while src[i].isspace():
                                i += 1
                                if i >= len(src):
                                    return ""
                            if src[i] == '{':
                                n = 1
                                i += 1
                                c = ''
                                while n != 0:
                                    c += src[i]
                                    i += 1
                                    if i >= len(src):
                                        return ""
                                    if src[i] == '{' and src[i - 1] != '\\':
                                        n += 1
                                    elif src[i] == '}' and src[i - 1] != '\\':
                                        n -= 1
                                args.append(latex2utf(c))
                            elif src[i] == '(':
                                n = 1
                                i += 1
                                c = '('
                                while n != 0:
                                    c += src[i]
                                    i += 1
                                    if i >= len(src):
                                        return ""
                                    if src[i] == '(' and src[i - 1] != '\\':
                                        n += 1
                                    elif src[i] == ')' and src[i - 1] != '\\':
                                        n -= 1
                                c += ')'
                                args.append(latex2utf(c))
                            else:
                                args.append(latex2utf(src[i]))
                        res += f.replacement(args)
                else:
                    res += f.replacement
                    
        elif src[i] in ['+', '-', '=', '<', '>']:
            if len(res) == 0 or not res[-1].isspace():
                res += ' '
            res += src[i] + ' '
        elif src[i] in ['^', '_'] and (i == 0 or src[i - 1] != '\\'):
            if src[i] == '^':
                f = cmd_mod
            else:
                f = cmd_sub
            i += 1
            while src[i].isspace():
                i += 1
                if i >= len(src):
                    return ""
            res = sub('\\s+$', '', res)
            if src[i] == '{':
                n = 1
                i += 1
                c = ''
                while n != 0:
                    c += src[i]
                    i += 1
                    if i >= len(src):
                        return ""
                    if src[i] == '{' and src[i - 1] != '\\':
                        n += 1
                    elif src[i] == '}' and src[i - 1] != '\\':
                        n -= 1
                res += f(latex2utf(c))
            elif src[i] == '(':
                n = 1
                i += 1
                c = '('
                while n != 0:
                    c += src[i]
                    i += 1
                    if i >= len(src):
                        return ""
                    if src[i] == '(' and src[i - 1] != '\\':
                        n += 1
                    elif src[i] == ')' and src[i - 1] != '\\':
                        n -= 1
                c += ')'
                res += f(latex2utf(c))
            else:
                res += f(latex2utf(src[i]))
        elif src[i] in ['{', '}' ] and not (i > 0 and src[i - 1] == '\\'):
            i += 1
            continue

        elif not src[i].isspace():
            res += src[i]

        i += 1

    res = sub('^\\s*', '', res)
    res = sub('\\s*$', '', res)
    res = sub('\\s\\s*', ' ', res)

    return res

if __name__ == '__main__':
    for i in stdin.readlines():
        print latex2utf(i)

# vim: set tabstop=4 softtabstop=4 shiftwidth=4 expandtab :

