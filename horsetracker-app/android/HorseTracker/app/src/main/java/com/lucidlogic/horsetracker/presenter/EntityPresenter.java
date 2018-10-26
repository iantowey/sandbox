package com.lucidlogic.horsetracker.presenter;

import com.lucidlogic.horsetracker.view.EntityView;

/**
 * Created by itowey on 11/08/16.
 */
public interface EntityPresenter {
    void setView(EntityView entityView);
    void updateView(String entityObjectId);
}
