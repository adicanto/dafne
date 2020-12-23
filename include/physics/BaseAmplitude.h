#pragma once

#include <string>

#include <hydra/detail/Config.h>

namespace dafne {

class BaseAmplitude
{
private:
	std::string _name;
	std::string _label;
	int _color;
	int _style;
	bool _isRemoved;

public:
	BaseAmplitude() = delete;

	__hydra_dual__ inline
	BaseAmplitude(const char *name="", const char *label="", int color=1, int style=1) : _name(name), _label(label), _color(color), _style(color), _isRemoved(false)
	{}

	__hydra_dual__ inline
	BaseAmplitude(BaseAmplitude const& other) : _name(other.Name()), _label(other.Label()), _color(other.Color()), _style(other.Style()), _isRemoved(other.IsRemoved())
	{}

	__hydra_dual__ inline
	BaseAmplitude& operator=(BaseAmplitude const& other)
	{
		if(this==&other) return *this;
		_name = other.Name();
		_label = other.Label();
		_color = other.Color();
		_style = other.Style();
		_isRemoved = other.IsRemoved();
		return *this;
	}

	__hydra_dual__ inline
	std::string Name() const { return _name; }

	__hydra_dual__ inline
	void SetName(const char *name) { _name = name; }

	__hydra_dual__ inline
	std::string Label() const { return _label; }
	
	__hydra_dual__ inline
	void SetLabel(const char *label) { _label = label; }

	__hydra_dual__ inline
	void Remove(bool value=true) { _isRemoved = value; }

	__hydra_dual__ inline
	bool IsRemoved() const { return _isRemoved; }

	__hydra_dual__ inline
	int Color() const { return _color; }

	__hydra_dual__ inline
	void SetColor(const int color) { _color = color; }

	__hydra_dual__ inline
	int Style() const { return _style; }
	
	__hydra_dual__ inline
	void SetStyle(const int style) { _style = style; }
};

} // namespace dafne
